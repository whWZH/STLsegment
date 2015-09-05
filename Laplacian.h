#include "STLModel.h"
#include <Eigen/Sparse>
typedef Eigen::SparseMatrix<double> SpMat; 
typedef Eigen::Triplet<double> T;
class CLaplacian
{
public:
	void change(vec_PVERT& m_vecLpvert,CSTLModel* &pCSTLModel,vector<POINT3D>& m_vecLp,vec_VECTOR3D& m_vecV);
	void back(vec_PVERT& m_vecLpvert,CSTLModel* &pCSTLModel,vector<POINT3D>& m_vecLp,vec_VECTOR3D& m_vecV);
	void selectP(CSTLModel* &pCSTLModel,vector<int>& m_vecSL);//Ñ¡Ôñ¿ØÖÆ¶¥µã
	void draw(COpenGLDC* pDC);
protected:
private:
};