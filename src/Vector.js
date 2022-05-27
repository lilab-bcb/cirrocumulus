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

  isSparse() {
    return false;
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
      },
    };
  }
}
