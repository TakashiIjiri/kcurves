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
	/// MainForm の概要
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
		bool isClosed()
		{
			return m_checkbox_closed->Checked;
    }

	protected:
		/// <summary>
		/// 使用中のリソースをすべてクリーンアップします。
		/// </summary>
		~MainForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::CheckBox^ m_checkbox_closed;
	private: System::Windows::Forms::CheckBox^ m_cb_catmullrom;
	protected:

	protected:

	protected:

	private:
		/// <summary>
		/// 必要なデザイナー変数です。
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// デザイナー サポートに必要なメソッドです。このメソッドの内容を
		/// コード エディターで変更しないでください。
		/// </summary>
		void InitializeComponent(void)
		{
			this->m_checkbox_closed = (gcnew System::Windows::Forms::CheckBox());
			this->m_cb_catmullrom = (gcnew System::Windows::Forms::CheckBox());
			this->SuspendLayout();
			// 
			// m_checkbox_closed
			// 
			this->m_checkbox_closed->AutoSize = true;
			this->m_checkbox_closed->Checked = true;
			this->m_checkbox_closed->CheckState = System::Windows::Forms::CheckState::Checked;
			this->m_checkbox_closed->Location = System::Drawing::Point(13, 13);
			this->m_checkbox_closed->Name = L"m_checkbox_closed";
			this->m_checkbox_closed->Size = System::Drawing::Size(59, 16);
			this->m_checkbox_closed->TabIndex = 0;
			this->m_checkbox_closed->Text = L"Closed";
			this->m_checkbox_closed->UseVisualStyleBackColor = true;
			this->m_checkbox_closed->CheckedChanged += gcnew System::EventHandler(this, &MainForm::m_checkbox_closed_CheckedChanged);
			// 
			// m_cb_catmulrom
			// 
			this->m_cb_catmullrom->AutoSize = true;
			this->m_cb_catmullrom->Checked = true;
			this->m_cb_catmullrom->CheckState = System::Windows::Forms::CheckState::Checked;
			this->m_cb_catmullrom->Location = System::Drawing::Point(13, 35);
			this->m_cb_catmullrom->Name = L"m_cb_catmulrom";
			this->m_cb_catmullrom->Size = System::Drawing::Size(154, 16);
			this->m_cb_catmullrom->TabIndex = 1;
			this->m_cb_catmullrom->Text = L"Show Catmul Rom α=0.5";
			this->m_cb_catmullrom->UseVisualStyleBackColor = true;
			// 
			// MainForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 12);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(684, 577);
			this->Controls->Add(this->m_cb_catmullrom);
			this->Controls->Add(this->m_checkbox_closed);
			this->DoubleBuffered = true;
			this->Margin = System::Windows::Forms::Padding(2);
			this->Name = L"MainForm";
			this->Text = L"MainForm";
			this->MouseDown += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::MainForm_MouseDown);
			this->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::MainForm_MouseMove);
			this->MouseUp += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::MainForm_MouseUp);
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion

	private: System::Void MainForm_MouseUp(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)  ;
	private: System::Void MainForm_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e);
	private: System::Void MainForm_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e);
	private: System::Void m_checkbox_closed_CheckedChanged(System::Object^ sender, System::EventArgs^ e);
	};


	inline void MainForm_repaint()
	{
		MainForm::getInst()->repaint();
	}
	inline bool MainForm_isClosed()
	{
    return MainForm::getInst()->isClosed();
	}

}
