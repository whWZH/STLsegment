#include "StdAfx.h"
#include "OpenGLDC.h"


COpenGLDC::COpenGLDC(HWND hWnd):m_hWnd(hWnd)
{
}


COpenGLDC::~COpenGLDC(void)
{
}

BOOL COpenGLDC::InitDC()
{
	if(m_hWnd == NULL)return FALSE;
	m_Camera.init();                   //初始照相机设置
	m_hDC =::GetDC(m_hWnd);            //获取设备场境

	PIXELFORMATDESCRIPTOR pfdWnd =
	{
		sizeof(PIXELFORMATDESCRIPTOR),//Structure size
		1,                            //Structure version number
		PFD_DRAW_TO_WINDOW|           //Property flags
		PFD_SUPPORT_OPENGL|           
		PFD_DOUBLEBUFFER,
		PFD_TYPE_RGBA,
		24,                           //24-bit color
		0,0,0,0,0,0,                  //Not concerned with these
		0,0,0,0,0,0,0,                //No alpha or accum buffer
		32,                           //32-bit depth buffer
		0,0,                          //No stencil or aux buffer
		PFD_MAIN_PLANE,               //Main layer type
		0,                            //Reserved
		0,0,0,                        //Unsupported

	};

	int pixelformat;
	if((pixelformat = ChoosePixelFormat(m_hDC,&pfdWnd))== 0)
	{
		AfxMessageBox("ChoosePixelFormat to wnd failed");
		return FALSE;
	}

	if (SetPixelFormat(m_hDC,pixelformat,&pfdWnd) == FALSE)
	{
		AfxMessageBox("SetPixelFormat failed");
	}

	m_hRC = wglCreateContext(m_hDC);           //创建渲染场境
	VERIFY(wglMakeCurrent(m_hDC,m_hRC));       //当前化渲染场境
	GLSetupRC();                               //初始化渲染场境
	wglMakeCurrent(NULL,NULL);                 //非当前化渲染场境
	return m_hRC!=0;

}

void COpenGLDC::GLResize(int w,int h)
{
	wglMakeCurrent(m_hDC,m_hRC);

	if(h==0)h=1;    //反正被0除
	if(w==0)w=1;
	m_Camera.set_screen(w,h);
	//wglMakeCurrent(NULL,NULL);
}

void COpenGLDC::GLSetupRC()
{
	m_bShading = TRUE;           //使用着色模式

	glEnable(GL_DEPTH_TEST);     //使用消影
	//glEnable(GL_CULL_FACE);      //不计算对象内部
    glFrontFace(GL_CCW);          //三角片顶点逆时针方向（CCW）旋转表示模型的外表面

	//设置光源
	GLfloat lightAmbient[] = {0.75f,0.75f,0.75f,1.0f}; //设置环境光的颜色组成
	GLfloat lightDiffuse[] = {1.0f,1.0f,1.0f,1.0f};    //漫反射光
	//GLfloat  specular[] = {1.0f,1.0f,1.0f,1.0f};     //镜面光
	//GLfloat  lightPos[] = {1.0f,1.0f,1.0f,0.0f};     //光源位置，沿矢量（1,1,1）方向无穷远处

	glEnable(GL_LIGHTING);                           //使用光照模式
	glLightfv(GL_LIGHT0,GL_AMBIENT,lightAmbient);    //为光源0设置环境光
	glLightfv(GL_LIGHT0,GL_DIFFUSE,lightDiffuse);    //为光源0设置漫反射光
	//glLightfv(GL_LIGHT0,GL_SPECULAR,specular);       //为光源0设置镜面光
	//glLightfv(GL_LIGHT0,GL_POSITION,lightPos);       //设置光源0位置
	SetLightDirection(1,1,1);                        //设置光源的方向
	glEnable(GL_LIGHT0);                             //打开光源0

	//设置默认的颜色属性
	//SetBkColor(RGB(0,0,0));            //设置背景的默认颜色（黑色）
	SetBkColor(RGB(240,255,255));            //设置背景的默认颜色（白色）
	SetMaterialColor(RGB(225,175,22)); //设置材料的默认颜色
	SetColor(RGB(255,0,255));        //设置框架显示的默认颜色（白色）
	SetHighlightColor(RGB(255,0,0));   //设置高亮度显示的颜色（红色）

	glPointSize(3.0);                 //设置点的绘制尺寸

	////设置材料性质
	//glColorMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE);
	//glMaterialfv(GL_FRONT,GL_SPECULAR,specular);
	//glMateriali(GL_FRONT,GL_SHININESS,100);

	////设置背景色，黑色
	//glClearColor(0.0f,0.0f,0.0f,1.0f);
	//glColor3ub(0,0,255);             //设置材料的默认颜色

}

void COpenGLDC::Ready()
{
	wglMakeCurrent(m_hDC,m_hRC);                      //当前化渲染场境
	ClearBKground();
	OnShading();
	m_Camera.projection();                            //绘制图形之前的取景
	//glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); //清除颜色缓存与深度缓存
}

void COpenGLDC::Finish()
{
	glFlush();         
	SwapBuffers(m_hDC);         //交换帧存
	wglMakeCurrent(m_hDC,NULL);  //非当前化渲染场境
}

/////////////////////////光照和材料的设置//////////////////////
void COpenGLDC::ClearBKground()                                 //背景刷新
{
	GLclampf r,g,b;
	//获取背景色变量m_clrBK的颜色RGB分量
	r = (GLclampf)GetRValue(m_clrBk)/255.0;
	g = (GLclampf)GetGValue(m_clrBk)/255.0;
	b = (GLclampf)GetBValue(m_clrBk)/255.0;

	glClearColor(r,g,b,0.0f);      //设置清屏的RGBA颜色
	//glClearColor(1.0f, 1.0f, 1.0f, 1.0f );
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); //清除颜色缓存与深度缓存
}

void COpenGLDC::OnShading()                                     //光照模式设定
{
	if (m_bShading)
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	} 
	else
	{
			glDisable(GL_LIGHTING);
			glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

			glEnable(GL_CULL_FACE);//开启剪裁模式
			glCullFace(GL_BACK);//消除背面
	}
}
void COpenGLDC::Shading(BOOL bShading)
{
	m_bShading = bShading;
}
BOOL COpenGLDC::IsShading()
{
	return m_bShading;
}

void COpenGLDC::Lighting(BOOL bLighting)                        //光照开关
{
	if(bLighting)
		glEnable(GL_LIGHTING);
	else
		glDisable(GL_LIGHTING);
}
BOOL COpenGLDC::IsLighting()
{
	GLboolean bLighting;
	glGetBooleanv(GL_LIGHTING,&bLighting);
	return bLighting;
}

void COpenGLDC::SetLightDirection(float dx,float dy,float dz) //设置光源方向
{
	m_vecLight[0] = dx;
	m_vecLight[1] = dy;
	m_vecLight[2] = dz;
	GLfloat lightPos[] = {dx,dy,dz,0.0f};
	glLightfv(GL_LIGHT0,GL_POSITION,lightPos);
}
void COpenGLDC::GetLightDirection(float& dx,float& dy,float& dz)
{
	dx = m_vecLight[0];
	dy = m_vecLight[1];
	dz = m_vecLight[2];
}

void COpenGLDC::SetMaterialColor(COLORREF clr)                   //设置材质
{
	m_clrMaterial = clr;    //将材料颜色保存在成员变量m_clrMaterial中
	BYTE r,g,b;
	r = GetRValue(clr);
	g = GetGValue(clr);
	b = GetBValue(clr);

	//将以上颜色设置为当前的材质属性
	GLfloat mat_amb_diff[] = {(GLfloat)r/255,(GLfloat)g/255,(GLfloat)b/255,1.0};
	glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,mat_amb_diff);
}
void COpenGLDC::GetMaterialColor(COLORREF& clr)
{
	clr = m_clrMaterial;
}

void COpenGLDC::SetBkColor(COLORREF clr)                  //设置背景颜色
{
	m_clrBk = clr;
}
void COpenGLDC::GetBkColor(COLORREF& clr)
{
	clr = m_clrBk;
}

void COpenGLDC::SetColor(COLORREF clr)                  //设置非光照模式下的绘图颜色
{
	m_clr = clr;
	BYTE r,g,b;
	r = GetRValue(clr);
	g = GetGValue(clr);
	b = GetBValue(clr);
	glColor3ub(r,g,b);  //设置绘图颜色
	//glColor3ub(0.0f,0.0f,0.0f);
}
void COpenGLDC::GetColor(COLORREF& clr)
{
	clr = m_clr;
}

void COpenGLDC::SetHighlightColor(COLORREF clr)         //设置高亮度颜色
{
	m_clrHighlight = clr;
}
void COpenGLDC::GetHighlightColor(COLORREF& clr)
{
	clr = m_clrHighlight;
}
void COpenGLDC::Highlight(BOOL bHighlight)
{
	BYTE r,g,b;
	if (bHighlight)
	{
		r = GetRValue(m_clrHighlight);
		g = GetGValue(m_clrHighlight);
		b = GetBValue(m_clrHighlight);
	} 
	else
	{
		r = GetRValue(m_clrMaterial);
		g = GetGValue(m_clrMaterial);
		b = GetBValue(m_clrMaterial);
	}

	GLfloat mat_amb_diff[] = {(GLfloat)r/255,(GLfloat)g/255,(GLfloat)b/255,1.0};
	glMaterialfv(GL_BACK,GL_AMBIENT_AND_DIFFUSE,mat_amb_diff);

}

///////////////////////////////////绘图函数/////////////////////////////////////////
void COpenGLDC::DrawPoint(const POINT3D& pt)   //绘制一个空间点，其大小在GLsetupRC定义
{
	glBegin(GL_POINTS);
	glPointSize(GLfloat(5));
	    glVertex3f(pt.x,pt.y,pt.z);
	glEnd();
}

void COpenGLDC::DrawLine(const POINT3D& sp,const POINT3D& ep)
{
	glBegin(GL_LINES);
	    glVertex3f(sp.x,sp.y,sp.z);
		glVertex3f(ep.x,ep.y,ep.z);
	glEnd();
}
void COpenGLDC::DrawLine2(const POINT3D& sp,const POINT3D& ep)
{
	glBegin(GL_LINES);
	glVertex3f(sp.x,sp.y,sp.z);
	glVertex3f(ep.x,ep.y,ep.z);
	glEnd();
}

void COpenGLDC::DrawTriChip(double n0,double n1,double n2,
	                     double v00,double v01,double v02,
						 double v10,double v11,double v12,
						 double v20,double v21,double v22)
{
	glBegin(GL_TRIANGLES);
	    glNormal3d(n0,n1,n2);
		glVertex3d(v00,v01,v02);
		glVertex3d(v10,v11,v12);
		glVertex3d(v20,v21,v22);
	glEnd();
}

void COpenGLDC::DrawCoord()
{
	BOOL bLighting = IsLighting();     //关闭光照模式
	Lighting(FALSE);

	double width,height;               //坐标轴的显示长度为视景体宽或高的20%
	m_Camera.get_view_rect(width,height);
	double len = min(width,height);
	len *= 0.2;

	POINT3D cPt,xPt,yPt,zPt;
	xPt.x = yPt.y =zPt.z = len;

	COLORREF old_clr;
	GetColor(old_clr);

	SetColor(RGB(255,0,0));           //X轴，红
	DrawLine(cPt,xPt);

	SetColor(RGB(0,255,0));           //Y轴，绿
	DrawLine(cPt,yPt);

	SetColor(RGB(0,0,255));           //Z轴，蓝
	DrawLine(cPt,zPt);

	Lighting(bLighting);            //恢复光照模式
	SetColor(old_clr);              //恢复绘图颜色

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

}
POINT3D COpenGLDC::BeginSelection(int xPos,int yPos)
{
	POINT3D P;
	m_bSlelectionMode=TRUE;
	wglMakeCurrent(m_hDC,m_hRC);
	GLint viewport[4];
	glSelectBuffer(BUFFER_LENGTH,m_selectBuff);
	P=m_Camera.selection(xPos,yPos);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	return P;
	//InitNames();
}
int COpenGLDC::EndSelection(UINT* items)
{
	m_bSlelectionMode=FALSE;
	int hits=glRenderMode(GL_RENDER);
	for (int i=0;i<hits;i++)
	{
		items[i]=m_selectBuff[i*4+3];
	}
	wglMakeCurrent(m_hDC,NULL);
	return hits;
}
BOOL COpenGLDC::isSelectionMode()
{
	return m_bSlelectionMode;
}
//void COpenGLDC::InitNames()
//{
//	glInitNames();
//}
//void COpenGLDC::LoadName(UINT name)
//{
//	glLoadName(name);
//}
////////////////////////////////////////////////////////////////

//void COpenGLDC::RenderInformation(int txtA,int txtB)
//{
//	glEnable(GL_COLOR_MATERIAL);//启用顶点颜色，没有启用的情况下glColor4f不起作用
//	glDisable(GL_LIGHTING);////作为背景色需要暂时关闭光照，避免其他物体散色光反射，从而影响背景矩形
//	glMatrixMode(GL_PROJECTION);
//	glPushMatrix();
//	glLoadIdentity();
//	glMatrixMode(GL_MODELVIEW);
//	glPushMatrix();
//	glLoadIdentity();
//	glColor3f(0,0,0);
//	glRasterPos2f(-0.99,-0.97);
//	DrawText("triangle:",12);
//	CString str;
//	str.Format("%d",txtA);
//	char *pBuff=str.GetBuffer(0);
//	DrawText(pBuff,12);
//	glRasterPos2f(-0.99,-0.90);
//	DrawText("vertice:",12);
//	CString vn;
//	vn.Format("%d",txtB);
//	pBuff=vn.GetBuffer(0);
//	DrawText(pBuff,12);
//	// 恢复原来的矩阵
//	glMatrixMode(GL_PROJECTION);
//	glPopMatrix();
//	glMatrixMode(GL_MODELVIEW);
//	glPopMatrix();
//	glDisable(GL_COLOR_MATERIAL);//启用材质自己的颜色，而非顶点颜色
//	glEnable(GL_LIGHTING);//打开光照
//}
////绘制文字
//void COpenGLDC::DrawText(char* string,int flag)
//{
//	char* p = NULL;
//
//	for (p = string; *p; p++)
//	{
//		if(flag==24)
//			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, *p);
//		else if(flag==18)
//			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *p);
//		else
//			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *p);
//	}
//}