cd objReader
mex -output objReaderMex -outdir ..\@TriangleMesh objReader.cpp tiny_obj_loader.cpp
cd ..\

cd NewtonAutodiffMex
mex -output NewtonMex -outdir ..\ computeHessianMex.cpp -I..\..\ -I..\..\libigl\include
cd ..\