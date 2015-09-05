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
	m_Camera.init();                   //��ʼ���������
	m_hDC =::GetDC(m_hWnd);            //��ȡ�豸����

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

	m_hRC = wglCreateContext(m_hDC);           //������Ⱦ����
	VERIFY(wglMakeCurrent(m_hDC,m_hRC));       //��ǰ����Ⱦ����
	GLSetupRC();                               //��ʼ����Ⱦ����
	wglMakeCurrent(NULL,NULL);                 //�ǵ�ǰ����Ⱦ����
	return m_hRC!=0;

}

void COpenGLDC::GLResize(int w,int h)
{
	wglMakeCurrent(m_hDC,m_hRC);

	if(h==0)h=1;    //������0��
	if(w==0)w=1;
	m_Camera.set_screen(w,h);
	//wglMakeCurrent(NULL,NULL);
}

void COpenGLDC::GLSetupRC()
{
	m_bShading = TRUE;           //ʹ����ɫģʽ

	glEnable(GL_DEPTH_TEST);     //ʹ����Ӱ
	//glEnable(GL_CULL_FACE);      //����������ڲ�
    glFrontFace(GL_CCW);          //����Ƭ������ʱ�뷽��CCW����ת��ʾģ�͵������

	//���ù�Դ
	GLfloat lightAmbient[] = {0.75f,0.75f,0.75f,1.0f}; //���û��������ɫ���
	GLfloat lightDiffuse[] = {1.0f,1.0f,1.0f,1.0f};    //�������
	//GLfloat  specular[] = {1.0f,1.0f,1.0f,1.0f};     //�����
	//GLfloat  lightPos[] = {1.0f,1.0f,1.0f,0.0f};     //��Դλ�ã���ʸ����1,1,1����������Զ��

	glEnable(GL_LIGHTING);                           //ʹ�ù���ģʽ
	glLightfv(GL_LIGHT0,GL_AMBIENT,lightAmbient);    //Ϊ��Դ0���û�����
	glLightfv(GL_LIGHT0,GL_DIFFUSE,lightDiffuse);    //Ϊ��Դ0�����������
	//glLightfv(GL_LIGHT0,GL_SPECULAR,specular);       //Ϊ��Դ0���þ����
	//glLightfv(GL_LIGHT0,GL_POSITION,lightPos);       //���ù�Դ0λ��
	SetLightDirection(1,1,1);                        //���ù�Դ�ķ���
	glEnable(GL_LIGHT0);                             //�򿪹�Դ0

	//����Ĭ�ϵ���ɫ����
	//SetBkColor(RGB(0,0,0));            //���ñ�����Ĭ����ɫ����ɫ��
	SetBkColor(RGB(240,255,255));            //���ñ�����Ĭ����ɫ����ɫ��
	SetMaterialColor(RGB(225,175,22)); //���ò��ϵ�Ĭ����ɫ
	SetColor(RGB(255,0,255));        //���ÿ����ʾ��Ĭ����ɫ����ɫ��
	SetHighlightColor(RGB(255,0,0));   //���ø�������ʾ����ɫ����ɫ��

	glPointSize(3.0);                 //���õ�Ļ��Ƴߴ�

	////���ò�������
	//glColorMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE);
	//glMaterialfv(GL_FRONT,GL_SPECULAR,specular);
	//glMateriali(GL_FRONT,GL_SHININESS,100);

	////���ñ���ɫ����ɫ
	//glClearColor(0.0f,0.0f,0.0f,1.0f);
	//glColor3ub(0,0,255);             //���ò��ϵ�Ĭ����ɫ

}

void COpenGLDC::Ready()
{
	wglMakeCurrent(m_hDC,m_hRC);                      //��ǰ����Ⱦ����
	ClearBKground();
	OnShading();
	m_Camera.projection();                            //����ͼ��֮ǰ��ȡ��
	//glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); //�����ɫ��������Ȼ���
}

void COpenGLDC::Finish()
{
	glFlush();         
	SwapBuffers(m_hDC);         //����֡��
	wglMakeCurrent(m_hDC,NULL);  //�ǵ�ǰ����Ⱦ����
}

/////////////////////////���պͲ��ϵ�����//////////////////////
void COpenGLDC::ClearBKground()                                 //����ˢ��
{
	GLclampf r,g,b;
	//��ȡ����ɫ����m_clrBK����ɫRGB����
	r = (GLclampf)GetRValue(m_clrBk)/255.0;
	g = (GLclampf)GetGValue(m_clrBk)/255.0;
	b = (GLclampf)GetBValue(m_clrBk)/255.0;

	glClearColor(r,g,b,0.0f);      //����������RGBA��ɫ
	//glClearColor(1.0f, 1.0f, 1.0f, 1.0f );
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); //�����ɫ��������Ȼ���
}

void COpenGLDC::OnShading()                                     //����ģʽ�趨
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

			glEnable(GL_CULL_FACE);//��������ģʽ
			glCullFace(GL_BACK);//��������
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

void COpenGLDC::Lighting(BOOL bLighting)                        //���տ���
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

void COpenGLDC::SetLightDirection(float dx,float dy,float dz) //���ù�Դ����
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

void COpenGLDC::SetMaterialColor(COLORREF clr)                   //���ò���
{
	m_clrMaterial = clr;    //��������ɫ�����ڳ�Ա����m_clrMaterial��
	BYTE r,g,b;
	r = GetRValue(clr);
	g = GetGValue(clr);
	b = GetBValue(clr);

	//��������ɫ����Ϊ��ǰ�Ĳ�������
	GLfloat mat_amb_diff[] = {(GLfloat)r/255,(GLfloat)g/255,(GLfloat)b/255,1.0};
	glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,mat_amb_diff);
}
void COpenGLDC::GetMaterialColor(COLORREF& clr)
{
	clr = m_clrMaterial;
}

void COpenGLDC::SetBkColor(COLORREF clr)                  //���ñ�����ɫ
{
	m_clrBk = clr;
}
void COpenGLDC::GetBkColor(COLORREF& clr)
{
	clr = m_clrBk;
}

void COpenGLDC::SetColor(COLORREF clr)                  //���÷ǹ���ģʽ�µĻ�ͼ��ɫ
{
	m_clr = clr;
	BYTE r,g,b;
	r = GetRValue(clr);
	g = GetGValue(clr);
	b = GetBValue(clr);
	glColor3ub(r,g,b);  //���û�ͼ��ɫ
	//glColor3ub(0.0f,0.0f,0.0f);
}
void COpenGLDC::GetColor(COLORREF& clr)
{
	clr = m_clr;
}

void COpenGLDC::SetHighlightColor(COLORREF clr)         //���ø�������ɫ
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

///////////////////////////////////��ͼ����/////////////////////////////////////////
void COpenGLDC::DrawPoint(const POINT3D& pt)   //����һ���ռ�㣬���С��GLsetupRC����
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
	BOOL bLighting = IsLighting();     //�رչ���ģʽ
	Lighting(FALSE);

	double width,height;               //���������ʾ����Ϊ�Ӿ�����ߵ�20%
	m_Camera.get_view_rect(width,height);
	double len = min(width,height);
	len *= 0.2;

	POINT3D cPt,xPt,yPt,zPt;
	xPt.x = yPt.y =zPt.z = len;

	COLORREF old_clr;
	GetColor(old_clr);

	SetColor(RGB(255,0,0));           //X�ᣬ��
	DrawLine(cPt,xPt);

	SetColor(RGB(0,255,0));           //Y�ᣬ��
	DrawLine(cPt,yPt);

	SetColor(RGB(0,0,255));           //Z�ᣬ��
	DrawLine(cPt,zPt);

	Lighting(bLighting);            //�ָ�����ģʽ
	SetColor(old_clr);              //�ָ���ͼ��ɫ

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
//	glEnable(GL_COLOR_MATERIAL);//���ö�����ɫ��û�����õ������glColor4f��������
//	glDisable(GL_LIGHTING);////��Ϊ����ɫ��Ҫ��ʱ�رչ��գ�������������ɢɫ�ⷴ�䣬�Ӷ�Ӱ�챳������
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
//	// �ָ�ԭ���ľ���
//	glMatrixMode(GL_PROJECTION);
//	glPopMatrix();
//	glMatrixMode(GL_MODELVIEW);
//	glPopMatrix();
//	glDisable(GL_COLOR_MATERIAL);//���ò����Լ�����ɫ�����Ƕ�����ɫ
//	glEnable(GL_LIGHTING);//�򿪹���
//}
////��������
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