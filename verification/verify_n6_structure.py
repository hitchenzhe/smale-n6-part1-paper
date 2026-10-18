"""Independent exact checks for the N=6 manuscript's combinatorial bookkeeping.

This verifies the stated partition of labels 1..117 and the complement-product
matrix used in Lemma products. It does not recreate Chang--Chen's diagram
classification or certify the local asymptotic arguments.
"""
from itertools import combinations

ALGORITHM_II = {
    1, 10, 21, 23, 24, 25, 26, 33, 34, 35, 37, 38, 39, 41, 42, 53, 54,
    56, 59, 60, 68, 70, 72, 83, 85, 87, 91, 92, 96, 97, 113,
}
REPORTED_EXCLUDED = {
    2, 3, 4, 20, 22, 31, 32, 36, 40, 44, 45, 49, 50, 51, 55, 57, 61,
    77, 78, 84, 86, 95, 99, 100, 105, 109, 111, 112,
}
HISTORICAL_EXCLUDED = {
    5, 6, 52, 64, 65, 67, 69, 71, 73, 74, 75, 79, 80, 81, 82, 89, 90,
}
POSITIVE_EMPTY = {116}
REPORTED_SURVIVORS = {
    7, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 27, 28, 29, 30, 43, 46,
    47, 48, 58, 62, 76, 101, 102, 103, 104, 106, 107, 108, 110, 114,
    115, 117,
}
DERIVED = {8, 63, 66, 88, 93, 94, 98}


def determinant_bareiss(matrix):
    """Exact determinant of an integer square matrix (Bareiss elimination)."""
    a = [list(map(int, row)) for row in matrix]
    n = len(a)
    sign, previous = 1, 1
    for k in range(n - 1):
        pivot_row = next((r for r in range(k, n) if a[r][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
        for i in range(k + 1, n):
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def rank_rational(matrix):
    """Exact rank using Fraction-free integer row elimination."""
    from fractions import Fraction

    a = [[Fraction(v) for v in row] for row in matrix]
    rows, cols = len(a), len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [v / scale for v in a[rank]]
        for r in range(rows):
            if r != rank and a[r][col]:
                factor = a[r][col]
                a[r] = [x - factor * y for x, y in zip(a[r], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def verify_partition():
    groups = [
        ALGORITHM_II,
        REPORTED_EXCLUDED,
        HISTORICAL_EXCLUDED,
        POSITIVE_EMPTY,
        REPORTED_SURVIVORS,
        DERIVED,
    ]
    names = ["Algorithm II", "reported exclusions", "historical exclusions",
             "positive-empty", "reported survivors", "derived"]
    for i, left in enumerate(groups):
        for j in range(i + 1, len(groups)):
            overlap = left & groups[j]
            assert not overlap, f"overlap {names[i]} / {names[j]}: {sorted(overlap)}"
    union = set().union(*groups)
    assert union == set(range(1, 118)), f"missing={sorted(set(range(1,118))-union)}; extra={sorted(union-set(range(1,118)))}"
    assert list(map(len, groups)) == [31, 28, 17, 1, 33, 7]
    assert len(REPORTED_SURVIVORS | DERIVED) == 40
    print("diagram labels: 117/117 covered; groups are disjoint (31, 28, 17, 1, 33, 7)")
    print("exceptional-mass obligations: 33 published + 7 derived = 40")


def verify_complement_matrix():
    pairs = list(combinations(range(1, 7), 2))
    matrix = [
        [int(p == q) + int(set(p).isdisjoint(q)) for q in pairs]
        for p in pairs
    ]
    assert len(pairs) == 15 and all(len(row) == 15 for row in matrix)
    assert all(matrix[i][j] == matrix[j][i] for i in range(15) for j in range(15))
    assert all(sum(row) == 7 for row in matrix)
    expected_nullities = {7: 1, -2: 5, 2: 9}
    for eigenvalue, multiplicity in expected_nullities.items():
        shifted = [
            [matrix[i][j] - (eigenvalue if i == j else 0) for j in range(15)]
            for i in range(15)
        ]
        assert 15 - rank_rational(shifted) == multiplicity
    determinant = determinant_bareiss(matrix)
    assert determinant == -114688
    print("complement-product matrix: 15x15 symmetric; row sum 7")
    print("eigenvalue multiplicities: 7:1, -2:5, 2:9")
    print(f"determinant: {determinant}")


if __name__ == "__main__":
    verify_partition()
    verify_complement_matrix()
