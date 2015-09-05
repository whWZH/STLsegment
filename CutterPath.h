#include "calcubase.h"
#include "solidPolygon.h"
#include "OpenGLDC.h"
#include "eigen/Eigen"
#include "OBB.h"
using namespace Eigen;
using Eigen::MatrixXd;
class CcutterPath :public CSolidPolygon
{
public:
	CcutterPath();
	~CcutterPath();
public:
	//////////////////��ʼ���߹켣
	void update_direct(OBB& theBOX, VECTOR3D &Tdirection);//����о����
	void get_plane(vector<PVERT>& Bloop,double& A,double& B,double& C);//��С�������
	void CUT_line(OBB& theBOX,PPOLYGON &firstPLN,vector<vec_VECTOR3D>& CCnomal,vec_PHEDGE &vec_PH,vector<vec_PFACETTRI>& RBOX);//ƽ���з���
	void cut_partilize(vector<vector<PFACETTRI>>& Rbox, vec_PFACETTRI& m_vecPFacetTri);//���·���

	inline void get_interPOT(PPOLYGON &firstPLN, PPOLYGON &nextPLN, vec_VECTOR3D &Cnomal, double &h, double &r, vec_PHEDGE &vec_PH, vec_PFACETTRI& facelist, VECTOR3D &Tdirection);//�õ���һȦƫ�õ�
	inline void dele_interPOT(PPOLYGON firstPLN, PPOLYGON &nextPLN, vec_PFACETTRI& facelist);//ɾ������
	inline void ccpath_snooth(CPlane& cutPlane, vec_PFACETTRI& facelist, PPOLYGON &nextPLN, vec_VECTOR3D &Cnomal, vec_PHEDGE &vec_PH);//�켣�Ż�
	void change(PFACETTRI& theFac);//�ı���Ƭת��
	void get_CCpath(OBB& theBOX, vec_PPOLYPOLYGON& m_vecPPolyPolygons, double h, double r, vector<vec_PFACETTRI>& RBOX);//���һ��Ҷ�ӷֿ�Ĺ켣



	//PFACETTRI  OffsetPcc(PFACETTRI pFacet, VECTOR3D StepDir, POINT3D CurrentP, POINT3D& OffsetP, double ScallopHeight, double CutterR);   // ��ƫ�õ�
	//POINT3D ProjectionP_Sample(PFACETTRI pFacet, VECTOR3D StepDir, POINT3D CurrentP, CPlane& plane, double ScallopHeight, double CutterR,  double Interval);
	//bool CcutterPath::PlanewithEdge(PHEDGE pHEdge, CPlane Plane, POINT3D& p);
};
