#pragma once
#include "afx.h"

#include "gl/GL.h"
#include "gl/GLU.h"
#include "gl/glut.h"
#include "gl/glaux.h"

#include "GCamera.h"
#include "GeomBase.h"
#include "mesh.h"

#define BUFFER_LENGTH 64

class COpenGLDC :public CObject
{
public:
	COpenGLDC(HWND hWnd);
	virtual~COpenGLDC(void);

private:
	HWND   m_hWnd;       //�������ڵľ��
	HGLRC  m_hRC;        //��Ⱦ�������
	HDC    m_hDC;        //�豸�������

	COLORREF m_clrBk;       //���ڱ���ɫ
	COLORREF m_clr;         //�ǹ���ģʽ�µ�ģ����ɫ
	COLORREF m_clrHighlight;//���ڸ�����ʾʱ��ģ����ɫ����ʰȡһ������ʱ��Ҫ��������ʾ
	COLORREF m_clrMaterial; //���ϵ���ɫ

	BOOL     m_bShading;    //�Ƿ������ɫ��ʾ
	GLfloat  m_vecLight[3]; //��Դ����

//    BOOL     m_bSelectionMode;            //��ǰ�Ƿ���ѡ��ģʽ
//	  GLuint   m_selectBuff[BUFFER_LENGTH]; //ѡ�񻺴���

public:
	GCamera  m_Camera;    //�����������ȡ������

protected:
	void ClearBKground(); //���������ɫ
	void OnShading();     //����/�ǹ���ģʽ����

public:
	BOOL InitDC();                    //��ʼ��
	void GLResize(int cx,int cy);     //��Ӧ���ڳߴ�仯
	void GLSetupRC();                 //������Ⱦ����

	void Ready();                     //����ǰ׼������
	void Finish();                    //������ͼ����

////////////////////�Թ�������ɫ�Ĳ�������//////////////////////

	void Shading(BOOL bShading);        //����/�ǹ���ģʽ�л�
	BOOL IsShading();                   //��ǰ�Ƿ�����ɫģʽ

	void Lighting(BOOL bLighting);      //�Ƿ�ʹ�ù�Դ
	BOOL IsLighting();
	BOOL m_bSlelectionMode;
	GLuint m_selectBuff[BUFFER_LENGTH];

	void SetLightDirection(float dx,float dy,float dz);  //�������ȡ��Դ����
	void GetLightDirection(float& dx,float& dy,float& dz);

	void SetMaterialColor(COLORREF clr);            //�������ȡ������ɫ
	void GetMaterialColor(COLORREF& clr);

	void SetBkColor(COLORREF rgb);                 //�������ȡ������ɫ
	void GetBkColor(COLORREF& rgb);

	void SetColor(COLORREF rgb);                  //�������ȡ�ǹ���ģʽ�µĻ�����ɫ
	void GetColor(COLORREF& rgb);

	void SetHighlightColor(COLORREF clr);         //�������ȡ��������ʾ����ɫ
	void GetHighlightColor(COLORREF& clr);

	void Highlight(BOOL bLight = TRUE);           //������/������ʾ�л�
	POINT3D BeginSelection(int xPos,int yPos);
	int EndSelection(UINT* items);
	BOOL isSelectionMode();
	/*void InitNames();
	void LoadName(UINT name);
	void PushName(UINT name);
	void PopName();
*/

////////////////////////����OpenGl����ͼԪ�ĸ߼���ͼ����//////////////////////////////
	void DrawPoint(const POINT3D&);            //����һ���ռ��
	void DrawCoord();                          //�����û�����ϵ
	void DrawLine(const POINT3D& sp,const POINT3D& ep);//����һ��ֱ��
	void DrawLine2(const POINT3D& sp,const POINT3D& ep);
	void DrawTriChip(double n0,double n1,double n2,double v00,double v01,double v02,       //����һ��������Ƭ
		             double v10,double v11,double v12,double v20,double v21,double v22);
	/////////////////////////////////
	/*static void RenderInformation(int txtA,int txtB);
	static  void DrawText(char* string, int flag=18);*/

};
