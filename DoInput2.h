#pragma once


// DoInput2 �Ի���

class DoInput2 : public CDialogEx
{
	DECLARE_DYNAMIC(DoInput2)

public:
	DoInput2(CWnd* pParent = NULL);   // ��׼���캯��
	virtual ~DoInput2();

// �Ի�������
	enum { IDD = IDD_DIALOG2 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

	DECLARE_MESSAGE_MAP()
public:
	double theANG;
};
