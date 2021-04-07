// R chisq.test

import {fisherTest} from './FisherExact';


it('test_1', () => {
    /*
   TeaTasting <-
matrix(c(3, 1, 1, 3),
       nrow = 2,
       dimnames = list(Guess = c("Milk", "Tea"),
                       Truth = c("Milk", "Tea")))
fisher.test(TeaTasting)
     */
    const result = fisherTest(3, 1, 1, 3);
    expect(result).toBeCloseTo(0.4857);
});


it('test_2', () => {
    /*
 > Convictions <- matrix(c(2, 10, 15, 3), nrow = 2,
    +                       dimnames =
        +                           list(c("Dizygotic", "Monozygotic"),
            +                                c("Convicted", "Not convicted")))
     */
    const result = fisherTest(2, 15, 10, 3);
    expect(result).toBeCloseTo(0.0005367);
});
