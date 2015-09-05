#pragma once
#include "..\trimesh\include\TriMesh.h"

class CToolFunction
{
public:
	CToolFunction(void);
	~CToolFunction(void);
    static bool NearestPointMeshV(vec vi,TriMesh*Tmesh,int &vIndex);
    static bool Extract_Beselect_Mesh(TriMesh * CTmesh,TriMesh*&ExtractMesh,vector<int>&Map_Extract_To_Base_V,vector<int>&Map_Extract_To_Base_F);
	
};
