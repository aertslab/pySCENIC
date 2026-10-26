from functools import partial
from multiprocessing import get_all_start_methods, get_context

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csc_matrix

from pyscenic.cli import arboreto_with_multiprocessing as grn


class TrackedArray(np.ndarray):
    reductions = 0

    def __reduce_ex__(self, protocol):
        type(self).reductions += 1
        return super().__reduce_ex__(protocol)


def worker_inputs(sparse=False):
    rng = np.random.RandomState(13)
    matrix = rng.normal(size=(32, 5))
    matrix[:, 3] = 2 * matrix[:, 0] + 0.1 * rng.normal(size=32)
    matrix[:, 4] = matrix[:, 1] - matrix[:, 2]
    if sparse:
        matrix = csc_matrix(matrix)
    names = ["TF1", "TF2", "TF3", "G1", "G2"]
    tf_matrix, tf_names = grn.to_tf_matrix(matrix, names, names[:3])
    return names, matrix, tf_matrix, tf_names


@pytest.mark.parametrize("start_method", get_all_start_methods())
@pytest.mark.parametrize("sparse", [False, True])
@pytest.mark.parametrize("method", ["RF", "GBM"])
def test_initialized_workers_preserve_seeded_inference(start_method, sparse, method):
    names, matrix, tf_matrix, tf_names = worker_inputs(sparse)
    original_matrix = matrix.copy()
    original_tf_matrix = tf_matrix.copy()
    kwargs = dict(grn.RF_KWARGS if method == "RF" else grn.SGBM_KWARGS)
    kwargs["n_estimators"] = 12
    args = (names, matrix, tf_matrix, tf_names, [method, kwargs], 42)
    # Repeat TF targets after other targets to detect retained state changes.
    targets = [0, 1, 2, 3, 4, 0, 2]
    expected = [grn.run_infer_partial_network(target, *args) for target in targets]

    with get_context(start_method).Pool(
        2, initializer=grn._init_worker, initargs=args
    ) as pool:
        actual = list(pool.imap(grn._run_infer_partial_network, targets, chunksize=1))

    for observed, baseline in zip(actual, expected):
        assert not observed.empty
        pd.testing.assert_frame_equal(observed, baseline, check_exact=True)
        assert not (observed["TF"] == observed["target"]).any()
    if sparse:
        assert (matrix != original_matrix).nnz == 0
        assert (tf_matrix != original_tf_matrix).nnz == 0
    else:
        np.testing.assert_array_equal(matrix, original_matrix)
        np.testing.assert_array_equal(tf_matrix, original_tf_matrix)


@pytest.mark.parametrize("start_method", get_all_start_methods())
def test_matrices_are_not_serialized_per_target(start_method):
    names, matrix, tf_matrix, tf_names = worker_inputs()
    args = (
        names,
        matrix.view(TrackedArray),
        tf_matrix.view(TrackedArray),
        tf_names,
        ["RF", dict(grn.RF_KWARGS, n_estimators=2)],
        42,
    )
    targets = list(range(len(names)))
    context = get_context(start_method)
    TrackedArray.reductions = 0
    with context.Pool(2) as pool:
        baseline = list(
            pool.imap(
                partial(
                    grn.run_infer_partial_network,
                    gene_names=args[0],
                    ex_matrix=args[1],
                    tf_matrix=args[2],
                    tf_matrix_gene_names=args[3],
                    method_params=args[4],
                    seed=args[5],
                ),
                targets,
                chunksize=1,
            )
        )
    assert TrackedArray.reductions == 2 * len(targets)

    TrackedArray.reductions = 0
    with context.Pool(2, initializer=grn._init_worker, initargs=args) as pool:
        actual = list(pool.imap(grn._run_infer_partial_network, targets, chunksize=1))
    # fork inherits inputs; spawn/forkserver serialize them once per worker.
    assert TrackedArray.reductions == (0 if start_method == "fork" else 4)
    for observed, expected in zip(actual, baseline):
        pd.testing.assert_frame_equal(observed, expected, check_exact=True)
