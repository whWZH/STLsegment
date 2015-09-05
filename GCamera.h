#pragma once

#include "GeomBase.h"

#define VIEW_FRONT         0
#define VIEW_BACK          1
#define VIEW_TOP           2
#define VIEW_BOTTOM        3
#define VIEW_RIGHT         4
#define VIEW_LEFT          5
#define VIEW_SW_ISOMETRIC  6
#define VIEW_SE_ISOMETRIC  7
#define VIEW_NE_ISOMETRIC  8
#define VIEW_NW_ISOMETRIC  9

//#define ZOOM_ALL           10
//#define ZOOM_IN            11
//#define ZOOM_OUT           12

class GCamera
{
protected:
	//�ӵ�λ�ú͹۲췽��
	POINT3D      m_eye;
	POINT3D      m_ref;
	VECTOR3D     m_vecUp;

	//�Ӿ������
	double       m_far,m_near;
	double       m_width,m_height;

	//�ӿڲ���
	double       m_screen[2];
public:
	GCamera(void);
	~GCamera(void);
	//��ʼ������
	void init();
	POINT3D selection(int xPos,int yPos);//ѡȡģʽ
	//ȡ������,���û�ͼģʽ�µ�ͶӰ�任
	void projection();
	void BeginGetMatrix();
	GLdouble ModelMatrix[16]; 
	GLdouble ProjMatrix[16];

	//ѡ���壬����ѡ��ģʽ�µ�ͶӰ�任
	//void selection(int xPos,int yPos);

	//��������
	void zoom(double scale);
	void zoom_all(double x0,double y0,double z0,double x1,double y1,double z1);
	GLdouble ox, oy, oz;

	//����ƽ��
	void move_view(double dpx,double dpy);

	//ѡ�������ͼ
	void set_view_type(int type);

	//�����ӿڳߴ�
	void set_screen(int x,int y);

	//�����ӵ�λ�úͷ���
	void set_eye(double eye_x,double eye_y,double eye_z);
	void set_ref(double ref_x,double ref_y,double ref_z);
	void set_vecUp(double up_dx,double up_dy,double up_dz);

	//�����Ӿ���
	void set_view_ret(double width,double height);
	void get_view_rect(double& width,double& height);
protected:
	//���������Ϸ����ʸ��
	void update_upVec();
};

