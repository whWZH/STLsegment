#include "../stdafx.h "
#include "ToolFunction.h"
#include "..\trimesh\include\TriMesh_algo.h"
CToolFunction::CToolFunction(void)
{
}

CToolFunction::~CToolFunction(void)
{
}
bool CToolFunction::NearestPointMeshV(vec vi,TriMesh*Tmesh,int &vIndex)
{
	if(Tmesh==NULL) 
	{
		return false;
	}
	Tmesh->need_use_KDTree();
	float	*v0=&Tmesh->vertices[0][0];
	const float *match = Tmesh->m_Tmesh_KD->closest_to_pt(vi);
	if (match)
	{
		vIndex=(match - v0) / 3;
		return true;
	}
	return false;
}
bool CToolFunction::Extract_Beselect_Mesh(TriMesh * CTmesh,TriMesh*&ExtractMesh,vector<int>&Map_Extract_To_Base_V,vector<int>&Map_Extract_To_Base_F)
{   
	TriMesh*Tmesh=CTmesh;
	Tmesh->need_neighbors();
	Tmesh->need_adjacentedges();
	Tmesh->need_adjacentfaces();
	TriMesh*NewMesh =new TriMesh();
	int vn=Tmesh->vertices.size();
	int fn=Tmesh->faces.size();
	//处理内外角点
	vector<int>beselectF;
	beselectF.reserve(fn);
	vector<bool>insertflagf(fn,false);

	for (int i=0;i<fn;i++)
	{
		if (Tmesh->faces[i].beSelect==true)
		{   
			for (int j=0;j<3;j++)
			{   
				int vid=Tmesh->faces[i][j];
				int afn=Tmesh->adjacentfaces[vid].size();
				vector<int>&af=Tmesh->adjacentfaces[vid];
				int countaf=0;//计算邻接三角面片中被选中的三角面片
				int noselectf;
				for (int k=0;k<afn;k++)
				{   
					TriMesh::Face f0=Tmesh->faces[af[k]];
					if (f0.beSelect==true)
					{
						countaf++;
					}
					else noselectf=af[k];
				}
				if (countaf==1)continue;//外角点
				if (insertflagf[i]==false)
				{
					beselectF.push_back(i);
					insertflagf[i]=true;
				}
				//内角点
				if(countaf==(afn-1))
				{        
					if (insertflagf[noselectf]==false)
					{
						beselectF.push_back(noselectf);
						insertflagf[noselectf]=true;
					}		
				}
			}
		}

	}
	vector<bool>Visitflag(vn,false);
	vector<int>Vnewid(vn);
	TriMesh::Face newface;
	int fnbeselect=beselectF.size();
	if (fnbeselect==0)
	{
		AfxMessageBox("警告：未发现被选中的区域，或选中的区域过小");
		return false;
	}
	NewMesh->faces.resize(fnbeselect);
	NewMesh->vertices.reserve(vn);
	Map_Extract_To_Base_V.reserve(vn);
	Map_Extract_To_Base_F.resize(fnbeselect);
	for (int i=0;i<fnbeselect;i++)
	{
		for (int j=0;j<3;j++)
		{   
			int vid=Tmesh->faces[beselectF[i]][j];
			if (Visitflag[vid]==false)
			{   
				Vnewid[vid]=NewMesh->vertices.size();
				NewMesh->vertices.push_back(Tmesh->vertices[vid]);
				Map_Extract_To_Base_V.push_back(vid);//映射
				Visitflag[vid]=true;
			}

			newface[j]=Vnewid[vid];		
		}
		NewMesh->faces[i]=newface;
		Map_Extract_To_Base_F[i]=beselectF[i];


	}

	NewMesh->Is_clouds=false;
	ExtractMesh=NewMesh;
	return true;

}

