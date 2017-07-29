//
// Copyright 2012-2013, Syoyo Fujita.
//
// Licensed under 2-clause BSD liecense.
//
#pragma once
//#ifndef _TINY_OBJ_LOADER_H
//#define _TINY_OBJ_LOADER_H

#include <string>
#include <vector>
#include <map>

namespace tinyobj {

typedef struct
{
    std::string name;

    float ambient[3];
    float diffuse[3];
    float specular[3];
    float transmittance[3];
    float emission[3];
    float shininess;
    float ior;                // index of refraction

    std::string ambient_texname;
    std::string diffuse_texname;
    std::string specular_texname;
    std::string normal_texname;
    std::map<std::string, std::string> unknown_parameter;
} material_t;

typedef struct
{
	// nave: mayby add different vectors for triangulated vertices, textures and normals
	//std::vector<float>				  positions_triangulated;
	//std::vector<float>				  normals_triangulated;
	//std::vector<float>				  texcoords_triangulated;
    std::vector<float>				  positions;
    std::vector<float>				  normals;
    std::vector<float>				  texcoords;
    std::vector< std::vector<int> >   face_vertex_ind;
    std::vector< std::vector<int> >   face_normals_ind;
    std::vector< std::vector<int> >   face_texcoords_ind;
    std::vector<unsigned int>   indices;
	bool isTriangleMesh;
	bool isQuadMesh;
} mesh_t;

typedef struct
{
    std::string  name;
    material_t   material;
    mesh_t       mesh;
} shape_t;

/// Loads .obj from a file.
/// 'shapes' will be filled with parsed shape data
/// The function returns error string.
/// Returns empty string when loading .obj success.
/// 'mtl_basepath' is optional, and used for base path for .mtl file.
std::string LoadObj( std::vector<shape_t>& shapes,   // [output]
					 const char* filename,
					 const char* mtl_basepath = NULL);

};

//#endif  // _TINY_OBJ_LOADER_H
