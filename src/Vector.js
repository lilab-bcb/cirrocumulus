export class Vector {
    constructor(name, values) {
        this.name = name;
        this.values = values;
    }

    getName() {
        return this.name;
    }

    size() {
        return this.values.length;
    }

    get(i) {
        return this.values[i];
    }

    asArray() {
        return this.values;
    }


}
