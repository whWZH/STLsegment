#ifndef  CCircusSelection_H_
#define  CCircusSelection_H_
#pragma once
#include "..\stdafx.h"
#include "..\trimesh\include\TriMesh.h"
#include "..\trimesh\include\Vec.h"
#include "..\opengl\OpenGLDC.h"
#include "..\ViewSpace\ShowShader.h"
using namespace ViewSpace;
#include <set>
#include <list>
using namespace std;
namespace Sele{
class CCircusSelection
{
public:
	CCircusSelection(TriMesh*Tmesh,COpenGLDC *pDC);
	CCircusSelection(void);
	~CCircusSelection(void);
	void PenEnd();
    void PushCurvePoint(vec2 mouseposition2D);
	void Render(ShaderModel shadermodle0=SMOOTH_SHADER);





public:
	TriMesh* m_Tmesh;
	COpenGLDC* m_pGLDC;
	vector<vec2> m_PenCurve;
	vector<int>m_PenCurve_3D;
	vector<int>m_All_PenCurve_3D;
private:
	class CNode
	{    
		typedef CNode* CNode_pointer;
	public:
		CNode(void){m_VisitFlag=Inactive;m_GeodesicDistance =std::numeric_limits<float>::max();};
		~CNode(void){};
		enum VisitFlag{Frozen,Active,Inactive};
		VisitFlag m_VisitFlag;
		float m_GeodesicDistance;
		int  m_VertexID;
		int  m_Pre_Edge;

		float& distance_from_source(){return m_GeodesicDistance;};
		virtual bool operator()(CNode_pointer const s1, CNode_pointer const s2) const
		{
			return s1->distance_from_source()!=s2->distance_from_source() ?
				s1->distance_from_source() < s2->distance_from_source() ://不等的时候根据距离排序
			s1->m_VertexID< s2->m_VertexID;//相等的时候需要根据ID号来
		};
		void clear()
		{
			m_GeodesicDistance = std::numeric_limits<float>::max();//无穷大数值
			m_VisitFlag=Inactive;
		}
	};
	typedef CNode* CNode_pointer;


	void CircusSelection();
};
}
#endif
