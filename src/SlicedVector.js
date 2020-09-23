export class SlicedVector {
    constructor(vector, indices) {
        this.vector = vector;
        this.indices = indices;
    }

    getName() {
        return this.vector.getName();
    }

    size() {
        return this.indices.length;
    }

    get(i) {
        return this.vector.get(this.indices[i]);
    }
}
