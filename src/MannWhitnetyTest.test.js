import {mannWhitney} from './MannWhitneyUTest';
import {rankdata} from './util';

it('testMannWhitneyUSimple', () => {
    /* Target values computed using R version 2.11.1
           * x <- c(19, 22, 16, 29, 24)
           * y <- c(20, 11, 17, 12)
           * wilcox.test(x, y, alternative = "two.sided", mu = 0, paired = FALSE, exact = FALSE, correct = FALSE)
           * W = 17, p-value = 0.08641
           */
    const x = [19, 22, 16, 29, 24];
    const y = [20, 11, 17, 12];
    const result = mannWhitney(x, y);
    expect(result.p).toBeCloseTo(0.08641);
    expect(result.statistic).toBeCloseTo(17);
});

it('testBigDataSet', () => {
    const d1 = [];
    const d2 = [];
    for (let i = 0; i < 1500; i++) {
        d1[i] = 2 * i;
        d2[i] = 2 * i + 1;
    }
    const result = mannWhitney(d1, d2);
    expect(result.statistic).toBeGreaterThan(0.1);
});


it('testTiesRank', () => {
    // compare with scipy.stats.rankdata
    let ranks = rankdata([3, 3, 2, 1, 1]);
    let expected = [4.5, 4.5, 3, 1.5, 1.5];
    for (let i = 0; i < ranks.length; i++) {
        expect(ranks[i]).toBe(expected[i]);
    }

    ranks = rankdata([1, 1, 3, 3, 2]);
    expected = [1.5, 1.5, 4.5, 4.5, 3];
    for (let i = 0; i < ranks.length; i++) {
        expect(ranks[i]).toBe(expected[i]);
    }
});

it('testTies', () => {
    const result = mannWhitney([4, 3, 2], [3, 2, 1]);
    expect(result.statistic).toBe(7);
});

it('testSame', () => {
    const result = mannWhitney([3, 2, 1], [3, 2, 1]);
    expect(result.statistic).toBe(4.5);
    expect(result.p).toBe(1);
});

it('testDifferent', () => {
    const result = mannWhitney([6, 5, 4], [3, 2, 1]);
    expect(result.statistic).toBe(9.0);
    expect(result.p).toBeCloseTo(0.049534613435626706);
});

