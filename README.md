# C model of an iterative algorithm for calculating the exponential and logarithm functions.

## About

[exp-log.c](exp-log-c) is a C model of an iterative algorithm for calculating exponentials and logarithms of fixed-point numbers [Ch. 8, 1]. The iterative part is implemented in carry-save representation. The model can be used for an experimental study and design of bit-equivalent hardware implementations through parameterized fixed-point precision and number of iterations [Ch. 5, 2].

The algorithm implemented here can be used to calculate elementary functions of arguments in very narrow ranges, such as after range reduction of floating-pointing inputs [Ch. 11, 1], [Ch. 5, 2]. Range reduction is not included in this model.

The underlying fixed-point format is parameterized through the macro `FRACT_BITS` which defines the number of bits in the fractional part of the fixed-point format used, up to the limit of the total word length of 64 bits.

The number of iterations to run the iterative algorithm for is an input to the `exp_log_iterative` function.

## References

[1] Jean-Michel Muller. [Elementary Function: Algorithms and Implementation](https://www.springer.com/gp/book/9781489979810#aboutBook). Birkhäuser Basel. 3rd edition. 2016.

[2] Mantas Mikaitis. [Arithmetic Accelerators for a Digital Neuromorphic Processor](https://www.research.manchester.ac.uk/portal/files/173360010/FULL_TEXT.PDF). PhD thesis. University of Manchester. July 2020.

## Licence

See [license.txt](license.txt) for licensing information.
