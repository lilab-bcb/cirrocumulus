import {SparseVector} from './SparseVector';

it('iteration', () => {
    const v = new SparseVector('', [0, 3, 4], [1, 2, 3], 10);
    let count = 0;
    let sum = 0;
    for (const value of v) {
        count++;
        sum += value;
    }
    expect(sum).toBe(6);
    expect(count).toBe(3);
});
it('randomAccess', () => {
    const v = new SparseVector('', new Uint8Array([0, 3, 4]), [1, 2, 3], 10);
    const expectedMap = {0: 1, 3: 2, 4: 3};
    for (let i = 0; i < v.size(); i++) { // ascending
        const expectedValue = expectedMap[i] || 0;
        expect(v.get(i)).toBe(expectedValue);
    }
    for (let i = v.size() - 1; i >= 0; i--) { // descending
        const expectedValue = expectedMap[i] || 0;
        expect(v.get(i)).toBe(expectedValue);
    }
});
