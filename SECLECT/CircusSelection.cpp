#include "..\stdafx.h"
#include "CircusSelection.h"
#include "OpenglSelection.h"
#include "..\BaseToolClass\dijkstrapath.h"
#include "..\BaseToolClass\ToolFunction.h"
using namespace Sele;
CCircusSelection::CCircusSelection(void)
{   

}
CCircusSelection::~CCircusSelection(void)
{
}
CCircusSelection::CCircusSelection(TriMesh*Tmesh,COpenGLDC *pDC)
{   
	m_Tmesh=Tmesh;
	m_pGLDC=pDC;
	m_Tmesh->need_across_edge();
	//m_Tmesh->need_adjacentedges();
	m_Tmesh->need_neighbors();
	m_Tmesh->need_normals();
	
}
//参数为屏幕点，与对应的三维空间点
void CCircusSelection::PushCurvePoint(vec2 mouseposition2D)
{
	m_PenCurve.push_back(mouseposition2D);
	vec mouseposition3D=COpenglSelection::GetMousePoint3D(m_pGLDC,CPoint(mouseposition2D[0],mouseposition2D[1]));
	int NearestV=-1;
	if (CToolFunction::NearestPointMeshV(mouseposition3D,m_Tmesh,NearestV))
	m_PenCurve_3D.push_back(NearestV);
}
void CCircusSelection::Render(ShaderModel shadermodle0)
{
	CShowShader::ShaderMesh(m_Tmesh,shadermodle0);
	//绘制画笔
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);
	glDisable(GL_DEPTH_TEST);//透明物体，不启用深度，这样鼠标拾取的时候可以避开透明物

	glDisable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);//启用顶点颜色，没有启用的情况下glColor4f不起作用
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();//(投影矩阵和模型视图矩阵在这里起什么作用，颠倒位置没关系？）

	glColor4f(255.0/255,255.0/255,255.0/255,0.45);//glColor4f(64.0/255,224.0/255,208.0/255,0.45)
	glLineWidth(3.0);
	glEnable (GL_LINE_SMOOTH);//开启直线反走样，默认情况下是关闭的
	glHint (GL_LINE_SMOOTH, GL_NICEST);//行为控制

    glColor3f(0,1,0);
	glBegin(GL_LINE_STRIP);
	int cn=m_PenCurve.size();
	for (int i=0;i<cn;i++)
	{
		float x0=float(m_PenCurve[i][0])/viewport[2];
		float y0=float(viewport[3]-m_PenCurve[i][1]- 1)/viewport[3];
		float xyscale=float(viewport[2])/viewport[3];
		//原点坐标转化到-1到1之间
		x0=2*x0-1;
		y0=2*y0-1;
		glVertex2f(x0,y0);
	}

	glEnd();
	glDisable(GL_LINE_SMOOTH);

	// 恢复原来的矩阵
	glLineWidth(1.0);
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glEnable(GL_DEPTH_TEST);//关闭深度测试功能，深度测试保证只绘制前面的，挡住的不绘制
	glEnable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);
	

}
//笔选，置于鼠标松开后
void CCircusSelection::PenEnd()
{  
	if(m_PenCurve_3D.size()<5)return;
	CircusSelection();//画笔圈选

	m_PenCurve.clear();
	m_PenCurve_3D.clear();
	m_All_PenCurve_3D.clear();
}
void CCircusSelection::CircusSelection()
{	
	vector<CNode>m_nodes(m_Tmesh->vertices.size());
	int vn=m_Tmesh->vertices.size();
	for(int i=0;i<vn;i++)
	{
		m_nodes[i].m_VertexID=i;
	}

	int cn=m_PenCurve_3D.size();
	vector<int> m_AllVetex;
 

	for (int i=0;i<cn-1;i++)
	{
		DijkstraPath m_path(m_Tmesh,m_PenCurve_3D[i],m_PenCurve_3D[i+1]);
		m_path.RunDijkstra();
		m_path.output(m_All_PenCurve_3D);
		int np=m_All_PenCurve_3D.size();
		for (int j=0;j<np;j++)
		{
			m_AllVetex.push_back(m_All_PenCurve_3D[j]);
		}
		m_All_PenCurve_3D.clear();
		
	}
	for (int i=0;i<cn;i++)
	{
		m_AllVetex.push_back(m_PenCurve_3D[i]);
	}
	
	int nPenCurve_3D=m_AllVetex.size();
	for (int i=0;i<nPenCurve_3D;i++)
	{
		vector<int>&a=m_Tmesh->adjacentfaces[m_AllVetex[i]];
		int n=a.size();
		for (int j=0;j<n;j++)
		{
			m_Tmesh->faces[a[j]].beSelect=true;
		}
	}

	//扩散

	vec2 average_point;//边界点的中心
	for (int i=0;i<cn;i++)
	{
		average_point=average_point+m_PenCurve[i];
	}
	average_point=float(1.0/cn)*average_point;
	vec average_Point3D=COpenglSelection::GetMousePoint3D(m_pGLDC,CPoint(average_point[0],average_point[1]));


	int Seed_Face=-1;
	int NearestV=-1;
	if (CToolFunction::NearestPointMeshV(average_Point3D,m_Tmesh,NearestV))
	{
		Seed_Face=m_Tmesh->adjacentfaces[NearestV][0];
	}
	else 
	{
		Seed_Face=0;
	}

	int count_selef=0;
	list<int>templist;
	templist.push_back(Seed_Face);
	int fn=m_Tmesh->faces.size();
	while(!templist.empty())
	{
		list<int>::iterator it=templist.begin();
		TriMesh::Face &across_edge=m_Tmesh->across_edge[*it];
		for (int j=0;j<3;j++)
		{
			int &fneighbor=across_edge[j];
			if(fneighbor<0||fneighbor>=fn)continue;
			if (!m_Tmesh->faces[fneighbor].beSelect)
			{
				templist.push_back(fneighbor);
				m_Tmesh->faces[fneighbor].beSelect=true;
				count_selef++;
			}
		}
		templist.erase(templist.begin());
	}

}

