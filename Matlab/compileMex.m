cd objReader
mex -output objReaderMex -outdir ..\@TriangleMesh objReader.cpp tiny_obj_loader.cpp
cd ..\

cd NewtonAutodiffMex
mex -output NewtonMex -outdir ..\ computeHessianMex.cpp -I..\..\ -IC:\Dropbox\Projects\AutoCutCpp\libigl\include -IC:\Dropbox\Projects\AutoCutCpp\libigl\external\nanogui\ext\eigen
cd ..\