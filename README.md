# cb_math

Bunch of mathematical structures generally used in simulations/games/graphics. 


Put VectorMath.cpp/.h and MatrixMath.cpp/.h in your code base to use.
Test code is in main.cpp.
Run make to compile and then run ./main.exe to run the test.

## FAQ
#### Why not further modularize this?
Because most of these classes are used in conjunction with others (eg. Matrix4 with Vector4) . Also helps to keep the code compact and with fewer files to compile.

#### Why no pass by reference for some of the Vector2/3/4 function parameters?
Performance gains. The compiler finds it hard to optimize around const references for these small-sized objects.

License: Public domain
