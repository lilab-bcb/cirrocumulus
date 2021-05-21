import {sortedIndexOf} from 'lodash';


export class SparseVector {
    constructor(name, indices, values, size) {
        this.name = name;
        if (indices.length !== values.length) {
            throw new Error();
        }
        this.indices = indices;
        this.values = values;
        this.n = size;
    }

    getName() {
        return this.name;
    }

    isSparse() {
        return true;
    }

    size() {
        return this.n;
    }

    get(i) {
        const index = sortedIndexOf(this.indices, i);
        return index === -1 ? 0 : this.values[index];
    }

    [Symbol.iterator]() {
        let index = 0;
        const size = this.values.length;
        return {
            next: () => {
                if (index < size) {
                    return {value: this.values[index++], done: false};
                } else {
                    return {done: true};
                }
            }
        };
    }
}
