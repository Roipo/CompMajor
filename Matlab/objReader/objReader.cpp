#include <iostream>
#include <mex.h>
#include "tiny_obj_loader.h"


void getV(mxArray* outputStruct, int VField, const tinyobj::shape_t& converter)
{
	int counter = 0, Nv = converter.mesh.positions.size() / 3;
	mxArray* VFieldValue = mxCreateDoubleMatrix(Nv,3,mxREAL);
	double* V = mxGetPr(VFieldValue);
	for(int v = 0; v < Nv; v++)
	{
		V[counter] = converter.mesh.positions[3*v + 0];
		counter++;
	}
	for(int v = 0; v < Nv; v++)
	{
		V[counter] = converter.mesh.positions[3*v + 1];
		counter++;
	}
	for(int v = 0; v < Nv; v++)
	{
		V[counter] = converter.mesh.positions[3*v + 2];
		counter++;
	}
    /* Use mxSetFieldByNumber instead of mxSetField for efficiency
     * mxSetField(plhs[0],i,"name",mxCreateString(friends[i].name); */
	mxSetFieldByNumber(outputStruct,0,VField,VFieldValue);
}

void getVt(mxArray* outputStruct, int VtField, const tinyobj::shape_t& converter)
{
	int counter = 0, uvSize = converter.mesh.texcoords.size() / 2;
	if(uvSize == 0)
	{
		mxSetFieldByNumber(outputStruct,0,VtField,0);
		return;
	}
	mxArray* VtFieldValue = mxCreateDoubleMatrix(uvSize,2,mxREAL);
	double* Vt = mxGetPr(VtFieldValue);
	for(int i = 0; i < uvSize; i++)
	{
		Vt[counter] = converter.mesh.texcoords[2*i + 0];
		counter++;
	}
	for(int i = 0; i < uvSize; i++)
	{
		Vt[counter] = converter.mesh.texcoords[2*i + 1];
		counter++;
	}
    /* Use mxSetFieldByNumber instead of mxSetField for efficiency
     * mxSetField(plhs[0],i,"name",mxCreateString(friends[i].name); */
	mxSetFieldByNumber(outputStruct,0,VtField,VtFieldValue);
}

void getVn(mxArray* outputStruct, int VnField, const tinyobj::shape_t& converter)
{
	int counter = 0, normalsSize = converter.mesh.normals.size() / 3;
	if(normalsSize == 0)
	{
		mxSetFieldByNumber(outputStruct,0,VnField,0);
		return;
	}
	mxArray* VnFieldValue = mxCreateDoubleMatrix(normalsSize,3,mxREAL);
	double* Vt = mxGetPr(VnFieldValue);
	for(int i = 0; i < normalsSize; i++)
	{
		Vt[counter] = converter.mesh.normals[3*i + 0];
		counter++;
	}
	for(int i = 0; i < normalsSize; i++)
	{
		Vt[counter] = converter.mesh.normals[3*i + 1];
		counter++;
	}
	for(int i = 0; i < normalsSize; i++)
	{
		Vt[counter] = converter.mesh.normals[3*i + 2];
		counter++;
	}
    /* Use mxSetFieldByNumber instead of mxSetField for efficiency
     * mxSetField(plhs[0],i,"name",mxCreateString(friends[i].name); */
	mxSetFieldByNumber(outputStruct,0,VnField,VnFieldValue);
}

void F2CellArray(mxArray* outputStruct, int FField, const tinyobj::shape_t& converter)
{
	int faceVerticesSize = converter.mesh.face_vertex_ind.size();
	mxArray* FFieldValue = mxCreateCellMatrix(1,(mwSize)faceVerticesSize);

	/* Fill cell matrix with input arguments */
    for( mwIndex i = 0; i < (mwIndex)faceVerticesSize; i++)
	{
		int size = converter.mesh.face_vertex_ind[i].size();
		mxArray* faceVerts = mxCreateDoubleMatrix(1,size,mxREAL);
		double* faceVertsPtr = mxGetPr(faceVerts);
		for(int j = 0; j < size; j++)
		{
			double temp = converter.mesh.face_vertex_ind[i][j];
			faceVertsPtr[j] = converter.mesh.face_vertex_ind[i][j];
		}
        mxSetCell(FFieldValue,i,faceVerts);
    }
	mxSetFieldByNumber(outputStruct,0,FField,FFieldValue);
}

void Ft2CellArray(mxArray* outputStruct, int FtField, const tinyobj::shape_t& converter)
{
	int faceTexturesSize = converter.mesh.face_texcoords_ind.size();
	mxArray* FtFieldValue = mxCreateCellMatrix(1,(mwSize)faceTexturesSize);

	/* Fill cell matrix with input arguments */
    for( mwIndex i = 0; i < (mwIndex)faceTexturesSize; i++)
	{
		int size = converter.mesh.face_texcoords_ind[i].size();
		mxArray* faceTextures = mxCreateDoubleMatrix(1,size,mxREAL);
		double* faceVertsPtr = mxGetPr(faceTextures);
		for(int j = 0; j < size; j++)
		{
			double temp = converter.mesh.face_texcoords_ind[i][j];
			faceVertsPtr[j] = converter.mesh.face_texcoords_ind[i][j];
		}
        mxSetCell(FtFieldValue,i,faceTextures);
    }
	mxSetFieldByNumber(outputStruct,0,FtField,FtFieldValue);
}

void Fn2CellArray(mxArray* outputStruct, int FnField, const tinyobj::shape_t& converter)
{
	int faceNoramlsSize = converter.mesh.face_normals_ind.size();
	mxArray* FnFieldValue = mxCreateCellMatrix(1,(mwSize)faceNoramlsSize);

	/* Fill cell matrix with input arguments */
    for( mwIndex i = 0; i < (mwIndex)faceNoramlsSize; i++)
	{
		int size = converter.mesh.face_normals_ind[i].size();
		mxArray* faceTextures = mxCreateDoubleMatrix(1,size,mxREAL);
		double* faceVertsPtr = mxGetPr(faceTextures);
		for(int j = 0; j < size; j++)
		{
			double temp = converter.mesh.face_normals_ind[i][j];
			faceVertsPtr[j] = converter.mesh.face_normals_ind[i][j];
		}
        mxSetCell(FnFieldValue,i,faceTextures);
    }
	mxSetFieldByNumber(outputStruct,0,FnField,FnFieldValue);
}

void getF(mxArray* outputStruct, int FField, const tinyobj::shape_t& converter)
{
	int counter = 0, faceVerticesSize = converter.mesh.face_vertex_ind.size();
	if(faceVerticesSize == 0)
	{
		mxSetFieldByNumber(outputStruct,0,FField,0);
		return;
	}
	int meshType;
	if(converter.mesh.isTriangleMesh)
	{
		meshType = 3;
	}
	else if(converter.mesh.isQuadMesh)
	{
		meshType = 4;
	}
	else
	{
		F2CellArray(outputStruct,FField,converter);
		return;
	}
	mxArray* FFieldValue = mxCreateDoubleMatrix(faceVerticesSize,meshType,mxREAL);
	double* F = mxGetPr(FFieldValue);
	for(int n = 0; n < converter.mesh.face_vertex_ind[0].size(); n++)
	{
		for(int i = 0; i < faceVerticesSize; i++)
		{
			F[counter] = converter.mesh.face_vertex_ind[i][n];
			counter++;
		}
	}
    /* Use mxSetFieldByNumber instead of mxSetField for efficiency
     * mxSetField(plhs[0],i,"name",mxCreateString(friends[i].name); */
	mxSetFieldByNumber(outputStruct,0,FField,FFieldValue);
}

void getFt(mxArray* outputStruct, int FtField, const tinyobj::shape_t& converter)
{
	int counter = 0, faceTexturesSize = converter.mesh.face_texcoords_ind.size();
	if(faceTexturesSize == 0)
	{
		mxSetFieldByNumber(outputStruct,0,FtField,0);
		return;
	}
	int meshType;
	if(converter.mesh.isTriangleMesh)
	{
		meshType = 3;
	}
	else if(converter.mesh.isQuadMesh)
	{
		meshType = 4;
	}
	else
	{
		Ft2CellArray(outputStruct,FtField,converter);
		return;
	}
	mxArray* FtFieldValue = mxCreateDoubleMatrix(faceTexturesSize,meshType,mxREAL);
	double* Ft = mxGetPr(FtFieldValue);
	for(int n = 0; n < converter.mesh.face_texcoords_ind[0].size(); n++)
	{
		for(int i = 0; i < faceTexturesSize; i++)
		{
			Ft[counter] = converter.mesh.face_texcoords_ind[i][n];
			counter++;
		}
	}
    /* Use mxSetFieldByNumber instead of mxSetField for efficiency
     * mxSetField(plhs[0],i,"name",mxCreateString(friends[i].name); */
	mxSetFieldByNumber(outputStruct,0,FtField,FtFieldValue);
}

void getFn(mxArray* outputStruct, int FnField, const tinyobj::shape_t& converter)
{
	int counter = 0, faceNormalsSize = converter.mesh.face_normals_ind.size();
	if(faceNormalsSize == 0)
	{
		mxSetFieldByNumber(outputStruct,0,FnField,0);
		return;
	}
	int meshType;
	if(converter.mesh.isTriangleMesh)
	{
		meshType = 3;
	}
	else if(converter.mesh.isQuadMesh)
	{
		meshType = 4;
	}
	else
	{
		Fn2CellArray(outputStruct,FnField,converter);
		return;
	}
	mxArray* FnFieldValue = mxCreateDoubleMatrix(faceNormalsSize,meshType,mxREAL);
	double* Fn = mxGetPr(FnFieldValue);
	for(int n = 0; n < converter.mesh.face_normals_ind[0].size(); n++)
	{
		for(int i = 0; i < faceNormalsSize; i++)
		{
			Fn[counter] = converter.mesh.face_normals_ind[i][n];
			counter++;
		}
	}
    /* Use mxSetFieldByNumber instead of mxSetField for efficiency
     * mxSetField(plhs[0],i,"name",mxCreateString(friends[i].name); */
	mxSetFieldByNumber(outputStruct,0,FnField,FnFieldValue);
}

#define NUMBER_OF_FIELDS 8

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
	const char *fieldNames[] = { "V" , "Vt" , "Vn" , "F" , "Ft" , "Fn" , "isTriangleMesh", "isQuadMesh"};

	char *inputBuf;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
		within an if statement, because it will never get to the else
		statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
		the MEX-file) */
	if(nrhs!=1)
	{
		mexErrMsgIdAndTxt( "MATLAB:invalidNumInputs",
				"One input required.");
	}
	else if(nlhs != 1) 
	{
		mexErrMsgIdAndTxt( "MATLAB:maxlhs",
				"One output required.");
	}

	/* input must be a string */
	if ( mxIsChar(prhs[0]) != 1)
	{
		mexErrMsgIdAndTxt( "MATLAB:inputNotString",
				"Input must be a string.");
	}

	/* input must be a row vector */
	if (mxGetM(prhs[0])!=1)
	{
		mexErrMsgIdAndTxt( "MATLAB:inputNotVector",
				"Input must be a row vector.");
	}
  
	/* copy the string data from prhs[0] into a C string input_ buf.    */
    inputBuf = mxArrayToString(prhs[0]);
	/* convert obj file to raw data*/
	std::vector<tinyobj::shape_t> shapes;
	std::string err = tinyobj::LoadObj(shapes, inputBuf);

	if (!err.empty()) {
		mexErrMsgIdAndTxt( "MATLAB:invalidFile",
				"Failed to open file: possible reasons\n1) obj file not in the same directory\n2) obj file needs mtl");
	}



	mxFree(inputBuf);

	/* Create a 1-by-1 array of structs. */ 
    plhs[0] = mxCreateStructMatrix(1, 1 , NUMBER_OF_FIELDS , fieldNames);
  
    /* This is redundant, but here for illustration.  Since we just
       created the structure and the field number indices are zero
       based, name_field will always be 0 and phone_field will always
       be 1 */
    int VField = mxGetFieldNumber(plhs[0],"V");
    int VtField = mxGetFieldNumber(plhs[0],"Vt");
	int VnField = mxGetFieldNumber(plhs[0],"Vn");
	int FField = mxGetFieldNumber(plhs[0],"F");
	int FtField = mxGetFieldNumber(plhs[0],"Ft");
	int FnField = mxGetFieldNumber(plhs[0],"Fn");
	int isTriangularField = mxGetFieldNumber(plhs[0],"isTriangleMesh");
	int isQuadField = mxGetFieldNumber(plhs[0],"isQuadMesh");

	getV(plhs[0],VField,shapes[0]);
	getVt(plhs[0],VtField,shapes[0]);
	getVn(plhs[0],VnField,shapes[0]);
	getF(plhs[0],FField,shapes[0]);
	getFt(plhs[0],FtField,shapes[0]);
	getFn(plhs[0],FnField,shapes[0]);

	mxArray* val1 = mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(val1) = shapes[0].mesh.isTriangleMesh;
	mxSetFieldByNumber(plhs[0],0,isTriangularField,val1);

	mxArray* val2 = mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(val2) = shapes[0].mesh.isQuadMesh;
	mxSetFieldByNumber(plhs[0],0,isQuadField,val2);
	return;
}
