#include "calcubase.h"
#include "OpenGLDC.h"
#include <queue>
#include "eigen/Eigen"
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include "fstream"
typedef Eigen::SparseMatrix<double> SpMat; 
typedef Eigen::Triplet<double> T;
using namespace Eigen;
using Eigen::MatrixXd;

class Cboundary
{
public:
	vec_PVERT vecPcut;
	Octree* pOctree;
	//////////////////////////////////////////////////////��ʽ��ȡ�㷨
	void featurPOINT(vec_PFACETTRI& m_vecPFacetTri,vec_PHEDGE& m_vecPHEdge,vec_PVERT& m_vecPVert,double &NVlimit);//�ҵ���������
	void featurPOINT_new(vec_PFACETTRI& m_vecPFacetTri,vec_PHEDGE& m_vecPHEdge,vec_PVERT& m_vecPVert);
	void Skeletonizing(vector<vector<PVERT>>& RBLOOP,vec_PFACETTRI& m_vecPFacetTri,vec_PHEDGE& m_vecPHEdge,vec_PVERT& m_vecPVert
		,vec_PEDGE& m_vecPEdge);//HU�ķ�����ȡ�Ǽ���
	void  Skeletonizing_new(vector<vector<PVERT>>& RBLOOP,vec_PFACETTRI& m_vecPFacetTri,vec_PVERT& m_vecPVert,vec_PHEDGE& m_vecPHEdge,vec_PEDGE& m_vecPEdge);//�����������ҹǼ���
	void  Skeletonizing_new_1(vector<vector<PVERT>>& RBLOOP,vec_PVERT& m_vecPVert);//������������
	void Skeletonizing_close(vector<vector<PVERT>>& RBLOOP,vec_PFACETTRI& m_vecPFacetTri,vec_PVERT& m_vecPVert,vec_PHEDGE& m_vecPHEdge,double &Ang);//�γɷ�յ�һȦ��
	void Skeletonizing_boundary(vec_PEDGE patchT,queue<PHEDGE>& patchB);//�ҵ�һƬ����ı߽��
	void Skeletonizing_link(vector<PVERT>& BLOOP,vec_PFACETTRI& m_vecPFacetTri,vec_PVERT& m_vecPVert,vec_PHEDGE& m_vecPHEdge,double &Ang);//ͨ���Ǽ��߱պ�
	void Skeletonizing_link_new(vector<PVERT>& BLOOP,vec_PFACETTRI& m_vecPFacetTri,vec_PVERT& m_vecPVert);//ͨ���Ǽ��߱պ�
	void Skeletonizing_smoothing(vector<vector<PFACETTRI>>& Rbox, vec_PFACETTRI& m_vecPFacetTri, vec_PVERT& m_vec, vector<vec_PVERT>& RBLOOP, vec_PHEDGE& m_vecph,
		vector<vector<bool>>& IDlist);//ȥ��ݹ�˳
	void change(PFACETTRI& theFac);//�ı���Ƭת��
	bool inverse(PFACETTRI theFac);//�ж�����true����ʱ��
	void texOUT(vec_PFACETTRI& upfacetri);//���STL�ļ�
	void resetLOOP(vector<vec_PVERT>& RBLOOP, vec_PVERT& m_vecPVert);//���¶���loop��
	void Skeletonizing_find_boundary(vec_PVERT ALOOP,queue<PHEDGE>& patchB,vec_PFACETTRI& m_vecPFacetTri,vec_PHEDGE& m_vecPHEdge,vec_PVERT& m_vecPVert);//�ҵ�һƬ����ı߽��
	BOOL Skeletonizing_delet(PHEDGE BIAN1,PHEDGE BIAN2);//�ж�Skeletonizing�߽��Ƿ���Ƴ�
	void Skeletonizing_partlize(vector<vector<PFACETTRI>>& Rbox,vec_PFACETTRI& m_vecPFacetTri,vec_PVERT& m_vecPVert);//������������������
	void Skeletonizing_partlize_new(vector<vector<PFACETTRI>>& Rbox, vec_PFACETTRI& m_vecPFacetTri);//������������������
	void Skeletonizing_combine(vector<vector<PFACETTRI>>& Rbox,vec_PFACETTRI& m_vecPFacetTri);
	//BOOL IntersectLinePlane(POINT3D LinePt,VECTOR3D LineNor, PFACETTRI pFacTri, POINT3D& InterP);//�����󽻵�
	void link_cone(vector<PVERT>& BLOOP,double angA,double angB,vec_PFACETTRI& m_vecPFacetTri,vec_PVERT& m_vecPVert,vec_PHEDGE& m_vecPHEdge);//����׶��������
	void sort_smooth(vec_PVERT& templine);
	/////////////////�����׾����ҹǼ���
	void Laboundary(vector<vector<PFACETTRI>>& Rbox,vector<vector<PVERT>>& RBLOOP,vec_PFACETTRI& m_vecPFacetTri,vec_PHEDGE& m_vecPHEdge,vec_PVERT& m_vecPVert);
	void Laboundary_L(vector<vector<PFACETTRI>>& Rbox,vec_PFACETTRI& facelist);
	void get_dis(PHEDGE A,double& theDis,double& theAng);//�õ������ͽǾ���
	void get_avg(double& avgDis,double& avgAng,vec_PFACETTRI& vec_face);//�õ�ƽ��ֵ
	void LaMeans(vector<vector<PFACETTRI>>& Rbox,vector<vector<double>>& theVector,vec_PFACETTRI& facelist);//K-MEANS����
	double LaMeans_dis(vector<double>& A,vector<double>& B);//����Kά��֮��ľ���
	/////////////////////////////////////////////������snak��˳
	void snke_parameter(vec_PVERT& vec_pot,vec_PFACETTRI& facelist);//������
	void snk_moving_2D(vector<POINT2D>& vec_paLine,vector<POINT2D>& vec_paAll);//��άsnake
	 /////////////////////////////////////
	 //////////////////////�Լ��ķ����ҹǼ���
	void snke_contract(vec_PVERT& Bloop,vec_PVERT& m_vecPVert);
	void snake_move(PVERT& Apot,PVERT& Ppot,PVERT& Npot);
	double snake_inside(PVERT& Apot,PVERT& Ppot,PVERT& Npot);
	double snake_exter(PVERT& Apot);
	BOOL snake_link(vector<PVERT>& Bloop,PVERT& ST,PVERT& EN,vector<PVERT>& m_vecPVert);
	void snake_smooth(vector<PVERT>& Bloop);//ȥ��ݹ�˳
	void snake_smooth_close(vector<PVERT>& Bloop);//ȥ��ݹ�˳
	void snake_delete(vector<vector<PVERT>>& deLoop);//
	double snak_energe(vector<PVERT>& Bloop);//����һ����������
	void snaking(vector<PVERT>& m_Pst,vector<PVERT>& m_vecPVert);
protected:
private:
};