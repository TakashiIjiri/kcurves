#pragma once

#include <vector>

namespace kcurves {


	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;





	/// <summary>
	/// MainForm �̊T�v
	/// </summary>
	public ref class MainForm : public System::Windows::Forms::Form
	{


		static MainForm^ m_singleton;

		MainForm(void)
		{
			InitializeComponent();
			
			System::Windows::Forms::Control::Paint += gcnew PaintEventHandler(this, &MainForm::RepaintFunction);
		}



	public:
		static MainForm^ getInst() {
			if (m_singleton == nullptr) {
				m_singleton = gcnew MainForm();
			}
			return m_singleton;
		}

		void repaint()
		{
			this->Refresh();
		}


		void RepaintFunction(Object^ sender, PaintEventArgs^ e);

	protected:
		/// <summary>
		/// �g�p���̃��\�[�X�����ׂăN���[���A�b�v���܂��B
		/// </summary>
		~MainForm()
		{
			if (components)
			{
				delete components;
			}
		}

	protected:

	protected:

	private:
		/// <summary>
		/// �K�v�ȃf�U�C�i�[�ϐ��ł��B
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// �f�U�C�i�[ �T�|�[�g�ɕK�v�ȃ��\�b�h�ł��B���̃��\�b�h�̓��e��
		/// �R�[�h �G�f�B�^�[�ŕύX���Ȃ��ł��������B
		/// </summary>
		void InitializeComponent(void)
		{
			this->SuspendLayout();
			// 
			// MainForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 12);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(426, 445);
			this->DoubleBuffered = true;
			this->Margin = System::Windows::Forms::Padding(2);
			this->Name = L"MainForm";
			this->Text = L"MainForm";
			this->MouseDown += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::MainForm_MouseDown);
			this->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::MainForm_MouseMove);
			this->MouseUp += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::MainForm_MouseUp);
			this->ResumeLayout(false);

		}
#pragma endregion

	private: System::Void MainForm_MouseUp(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)  ;
	private: System::Void MainForm_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e);
	private: System::Void MainForm_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e);
	};


	inline void MainForm_repaint()
	{
		MainForm::getInst()->repaint();
	}


}
