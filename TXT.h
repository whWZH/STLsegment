#include "calcubase.h"
#include "STLModel.h"
#include "SDF.h"
#include "gl\GL.h"
#include "gl\glaux.h"
#include "gl\GLU.h"
#include "gl\glut.h"
#include "GLView.h"
#include "KD-tree/KDtree.h"
#include "KD-tree/Vec.h"
#include <cstdlib>
#include <iostream>
class CTXT
{
public:
	CTXT(void);
	~CTXT(void);
   void	outputTXT(CSTLModel* pSTLModel);
   CSDF* pCSDF;
   void inputTXT(CSTLModel* pSTLModel);
protected:
private:
};