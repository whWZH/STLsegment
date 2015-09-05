#include "StdAfx.h"
#include <gl/GL.h>
#include <gl/GLU.h>

#include "GCamera.h"


GCamera::GCamera(void)
{
}
GCamera::~GCamera(void)
{
}
void GCamera::init()
{
	m_eye = POINT3D(0,0,1000.0);      //�ӵ�λ�ã�0��0,1000��
	m_ref = POINT3D(0,0,0);           //�ο���λ�ã�0��0��0��
	m_vecUp = VECTOR3D(0,1,0);         //�ӵ����Ϸ���0��1��0�� 

	m_far = 10000;
	m_near = 1;
	m_width = 2400.0;
	m_height = 2400.0;

	m_screen[0] = 400;
	m_screen[1] = 400;

}

void GCamera::projection()
{
	glMatrixMode(GL_PROJECTION);       //�л���ͶӰ�任��������
	glLoadIdentity();                  //��ʼ��ͶӰ����
	glRenderMode(GL_RENDER);           //��Ⱦģʽ

	//����ͶӰ����
	double left    =  -m_width/2.0;
	double right   =  m_width/2.0;
	double bottom  =  -m_height/2.0;
	double top     =  m_height/2.0;

	glOrtho(left,right,bottom,top,m_near,m_far);    //ֱ��ͶӰ��ʽ

	glMatrixMode(GL_MODELVIEW);                     //�л�����ͼ�任��������
	glLoadIdentity();                               //��ʼ����ͼ�任����

	//�����ӵ�λ�ú͹۲췽��
	gluLookAt(m_eye.x,m_eye.y,m_eye.z,m_ref.x,m_ref.y,m_ref.z,m_vecUp.dx,m_vecUp.dy,m_vecUp.dz);
	

}

void GCamera::set_screen(int x,int y)
{
	glViewport(0,0,x,y);
	if(y==0)y=1;
	double ratio = (double)x/(double)y;
	m_width *= (double)x/m_screen[0];
	m_height *= (double)y/m_screen[1];
	m_width = m_height*ratio;

	m_screen[0] = x;
	m_screen[1] = y;
}

void GCamera::set_eye(double eye_x,double eye_y,double eye_z)
{
	m_eye.x = eye_x;
	m_eye.y = eye_y;
	m_eye.z = eye_z;
}

void GCamera::set_ref(double ref_x,double ref_y,double ref_z)
{
	m_ref.x = ref_x;
	m_ref.y = ref_y;
	m_ref.z = ref_z;
}

void GCamera::set_vecUp(double up_dx,double up_dy,double up_dz)
{
	m_vecUp.dx = up_dx;
	m_vecUp.dy = up_dy;
	m_vecUp.dz = up_dz;
}

void GCamera::set_view_ret(double width,double height)
{
	m_width = width;
	m_height = height;
	double aspect = m_screen[0]/m_screen[1];
	m_width = m_height*aspect;
}

void GCamera::get_view_rect(double& width,double& height)
{
	width = m_width;
	height = m_height;
}

void GCamera::zoom(double scale)
{
	ASSERT(scale > 0.0);   //���ųߴ�������0
	m_height *= scale;     //�����Ӿ���ĸ�
	m_width  *= scale;     //�����Ӿ���Ŀ�
}

void GCamera::zoom_all(double x0,double y0,double z0,double x1,double y1,double z1)
{
	double width,height;
	double xx,yy,zz;

	//ģ�Ͱ��ݺеĳ����
	xx = x1-x0;
	yy = y1-y0;
	zz = z1-z0;

	//�����ܹ�����ģ�Ͱ��ݺе��Ӿ���Ŀ�͸�
	width = max(max(xx,yy),zz);
	height = max(max(xx,yy),zz);

	//���������Ӿ���Ŀ�͸�
	set_view_ret(width,height);

	//�ƶ��ӵ�Ͳο���
	VECTOR3D vec = m_eye - m_ref;
	m_ref.x = (x0+x1)/2;
	m_ref.y = (y0+y1)/2;
	m_ref.z = (z0+z1)/2;
	m_eye = m_ref + vec;
}

void GCamera::set_view_type(int type)   //���͹۲���ͼѡ��
{
	double r;
	VECTOR3D vec;

	vec = m_ref - m_eye;    //ʸ��vec��ʾ���߷���
	r = vec.GetLength();    //�ӵ�����յ�ľ���

	if(IS_ZERO(r)) r = 50.0;//��ֹ�ӵ�����յ��غ�
	if(r > 10000)  r = 10000;//��ֹ�ӵ������յ�̫Զ

	switch(type){
	case VIEW_FRONT:                       //ǰ��ͼ
		m_eye = m_ref+VECTOR3D(0,-r,0);    //�ƶ��ӵ�λ��
		m_vecUp = VECTOR3D(0,0,1);
		break;
	case VIEW_BACK:                       //����ͼ
		m_eye = m_ref+VECTOR3D(0,r,0);    //�ƶ��ӵ�λ��
		m_vecUp = VECTOR3D(0,0,1);
		break;
	case VIEW_TOP:                        //����ͼ
		m_eye = m_ref+VECTOR3D(0,0,r);    //�ƶ��ӵ�λ��
		m_vecUp = VECTOR3D(0,1,0);
		break;
	case VIEW_BOTTOM:                      //����ͼ
		m_eye = m_ref+VECTOR3D(0,0,-r);    //�ƶ��ӵ�λ��
		m_vecUp = VECTOR3D(0,1,0);
		break;
	case VIEW_RIGHT:                       //����ͼ
		m_eye = m_ref+VECTOR3D(r,0,0);    //�ƶ��ӵ�λ��
		m_vecUp = VECTOR3D(0,0,1);
		break;
	case VIEW_LEFT:                         //����ͼ
		m_eye = m_ref+VECTOR3D(-r,0,0);    //�ƶ��ӵ�λ��
		m_vecUp = VECTOR3D(0,0,1);
		break;
	case VIEW_SW_ISOMETRIC:                //SW ���ͼ
		m_eye = m_ref+VECTOR3D(-1,-1,1).GetNormal()*r;
		update_upVec();
		break;
	case VIEW_SE_ISOMETRIC:                //SE ���ͼ
		m_eye = m_ref+VECTOR3D(1,-1,1).GetNormal()*r;
		update_upVec();
		break;
	case VIEW_NE_ISOMETRIC:                //NE ���ͼ
		m_eye = m_ref+VECTOR3D(1,1,1).GetNormal()*r;
		update_upVec();
		break;
	case VIEW_NW_ISOMETRIC:                //NW ���ͼ
		m_eye = m_ref+VECTOR3D(-1,1,1).GetNormal()*r;
		update_upVec();
		break;
	}
}

void GCamera::update_upVec()
{
	VECTOR3D vec = m_ref - m_eye;  //���߷���ʸ��
	VECTOR3D zVec(0,0,1);
	VECTOR3D vec0;

	vec.Normalize();       //ʸ����λ��
	vec0 = vec*zVec;
	m_vecUp = vec0*vec;    //ʸ��m_vecUp�����߷���ֱ
}

void GCamera::move_view(double dpx,double dpy)
{
	VECTOR3D vec;
	VECTOR3D xUp,yUp;

	vec = m_ref - m_eye;        //���߷���ʸ��
	vec.Normalize();            //��λ��
	xUp = vec*m_vecUp;          //xUp�����ﴰ�ڵ�x���Ӧ��OpenGL�û�����ϵ��ʸ��
	yUp = xUp*vec;              //yUp�����ﴰ�ڵ�y���Ӧ��OpenGL�û�����ϵ��ʸ��

	m_eye -= xUp*m_width*dpx + yUp*m_height*dpy;//�ƶ��ӵ�λ��
	m_ref -= xUp*m_width*dpx + yUp*m_height*dpy;//�ƶ����յ�λ��
}
POINT3D GCamera::selection(int xPos,int yPos)
{
	POINT3D P;
	GLint vp[4];
	
	glGetIntegerv(GL_VIEWPORT,vp);//��ȡ��ǰ�ӿ���Ϣ
	//����ͶӰ�任
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//�л���ѡ��ģʽ
	glRenderMode(GL_SELECT);
	//�������λ�ã�����ѡ�����
	gluPickMatrix(xPos,vp[3]-yPos, 1, 1, vp );
	//����ͶӰ����
	double left		=  - m_width/2.0;
	double right	=  m_width/2.0;
	double bottom	=  - m_height/2.0;
	double top		=  m_height/2.0;
	 
	ox=oy=oz=0;
	glOrtho(left,right,bottom,top,m_near,m_far);
	GLint Viewport[4];
	glGetIntegerv(GL_VIEWPORT, Viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, ModelMatrix);		//���ģ�ͱ任����
	glGetDoublev(GL_PROJECTION_MATRIX, ProjMatrix);		//���ͶӰ�任����
	float fZValue = 0;
	glReadPixels(xPos,Viewport[3]-yPos,1,1,	GL_DEPTH_COMPONENT,GL_FLOAT,&fZValue );	
	
	gluUnProject(xPos, Viewport[3]-yPos,fZValue, ModelMatrix,ProjMatrix, Viewport,	&ox, &oy, &oz);
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity( );
	gluLookAt(m_eye.x,m_eye.y,m_eye.z,m_ref.x,m_ref.y,m_ref.z, m_vecUp.dx, m_vecUp.dy, m_vecUp.dz);
	P.x=ox;
	P.y=oy;
	P.z=oz;
	return P;
}

