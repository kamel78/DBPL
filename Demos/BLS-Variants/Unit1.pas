unit Unit1;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.ComCtrls, TreeView1, Unit2;

type
  TForm1 = class(TForm)
    TreeView11: TTreeView1;
    procedure TreeView11CustomDrawItem(Sender: TCustomTreeView; Node: TTreeNode;
      State: TCustomDrawState; var DefaultDraw: Boolean);
  private
    { Déclarations privées }
  public
    { Déclarations publiques }
  end;

var
  Form1: TForm1;

implementation

{$R *.dfm}

procedure TForm1.TreeView11CustomDrawItem(Sender: TCustomTreeView;
  Node: TTreeNode; State: TCustomDrawState; var DefaultDraw: Boolean);
begin
  with Sender as TCustomTreeView do
  begin
    if node.Level=0 then begin
                         Canvas.Font.Color := clBlue;
                         Canvas.Font.Style:=Canvas.Font.Style+[fsbold];
                         end;
      if node.Text[1]=''  then Canvas.Font.Name:='Symbol';
      if Node.ImageIndex=2 then begin
                                Canvas.Font.Color:=clred;
                                Canvas.Font.Style:=Canvas.Font.Style+[fsbold];
                                end;


  end;

end;

end.
