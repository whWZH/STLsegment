#pragma once


// DoInput �Ի���

class DoInput : public CDialogEx
{
	DECLARE_DYNAMIC(DoInput)

public:
	DoInput(CWnd* pParent = NULL);   // ��׼���캯��
	virtual ~DoInput();

// �Ի�������
	enum { IDD = IDD_DIALOG1 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

	DECLARE_MESSAGE_MAP()
public:
	double theLIMIT;
};
