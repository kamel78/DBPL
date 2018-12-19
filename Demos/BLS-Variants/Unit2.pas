unit Unit2;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.ComCtrls,  GeneralTypes,
  BLS24Curves, Vcl.StdCtrls,  Vcl.ExtCtrls, System.ImageList,System.IOUtils,EncdDecd,
  Vcl.ImgList, LargeIntegers, BNCurves, KSS18Curves,HashFunctions, ECCFp2,ECCFp3,ECCFp4,ECCFp,
  Vcl.Buttons;

type
  TForm2 = class(TForm)
    BLS24CurvePairing1: TBLS24CurvePairing;
    Panel1: TPanel;
    Panel2: TPanel;
    Edit1: TEdit;
    Edit2: TEdit;
    Edit3: TEdit;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Panel3: TPanel;
    Memo1: TMemo;
    Label9: TLabel;
    Label10: TLabel;
    Edit4: TEdit;
    Label11: TLabel;
    Image1: TImage;
    Button2: TButton;
    Button3: TButton;
    imgToolbar: TImageList;
    BNCurvePairing1: TBNCurvePairing;
    KSS18CurvePairing1: TKSS18CurvePairing;
    Button7: TButton;
    RichEdit1: TRichEdit;
    Button1: TButton;
    Button4: TButton;
    Button5: TButton;
    Panel4: TPanel;
    Label13: TLabel;
    Label12: TLabel;
    ComboBox1: TComboBox;
    ComboBox2: TComboBox;
    SaveDialog1: TSaveDialog;
    OpenDialog1: TOpenDialog;
    SpeedButton1: TSpeedButton;
    SpeedButton2: TSpeedButton;
    procedure WndProc(var Msgs: TMessage); override;
    procedure FormCreate(Sender: TObject);
    procedure Split(const Delimiter: Char; Input: string; const Strings: TStrings) ;
    procedure Button1Click(Sender: TObject);
    procedure Button6Click(Sender: TObject);
    Procedure BLS1Generate;
    Procedure BLS2Generate;
    Procedure BLS3Generate;
    Procedure BLS1Sign;
    Procedure BLS2Sign;
    Procedure BLS3Sign;
    Procedure BLS1Verify;
    Procedure BLS2Verify;
    Procedure BLS3Verify;
    Procedure EnableInterface;
    procedure Button2Click(Sender: TObject);
    procedure Button3Click(Sender: TObject);
    procedure Button4Click(Sender: TObject);
    procedure Button5Click(Sender: TObject);
    procedure SpeedButton2Click(Sender: TObject);
    procedure SpeedButton1Click(Sender: TObject);
  private
    { Déclarations privées }
  public
    { Déclarations publiques }
  end;
Type BLSwithoutROSecretKey=record
                        x,y:LargeInt;
                        end;
     BLSwithoutROPublicKey=record
                        case Seclevel:byte of
                        0:(u_Fp2,v_Fp2:G2BN);
                        1:(u_Fp3,v_Fp3:G2KSS18);
                        2:(u_Fp4,v_Fp4:G2BLS24);
                        end;
var
  Form2: TForm2;

  DefaultG2_128BN:G2BN;
  DefaultG1_128BN:G1BN;
  DefaultG2_192KSS:G2KSS18;
  DefaultG1_192KSS:G1KSS18;
  DefaultG2_256BLS:G2BLS24;
  DefaultG1_256BLS:G1BLS24;

  PublicKey_128BN:G2BN;
  PublicKey_192KSS:G2KSS18;
  PublicKey_256NBLS:G2BLS24;
  SecretKey:LargeInt;
  SecretKey_BLSWithoutRO:BLSwithoutROSecretKey;
  PublicKey_BLSWithoutRO:BLSwithoutROPublicKey;

  GeneratorsPairingsForBLS2_BN:GTBN;
  GeneratorsPairingsForBLS2_KSS:GTKSS18;
  GeneratorsPairingsForBLS2_BLS24:GtBls24;

  Tim:TTime;



implementation
 uses shellapi, richedit, Unit1, Unit3, Unit4;

{$R *.dfm}

function TimetoMessage(t:TTime):string;
var h,m,s,ms:word;
begin
DecodeTime(t,h,m,s,ms);
Result:=inttostr(h)+':'+inttostr(m)+':'+inttostr(s)+':'+inttostr(ms)+'ms';
end;

procedure TForm2.Button1Click(Sender: TObject);
begin
Form3.ShowModal;
if Form3.ModalResult=1 then begin
                            RichEdit1.Text:=Form3.RichEdit1.Text;
                            ComboBox1.ItemIndex:=Form3.ComboBox1.ItemIndex;
                            ComboBox2.ItemIndex:=Form3.ComboBox2.ItemIndex;
                            case Form3.ComboBox2.ItemIndex of
                            0:Label4.Caption:='Public Key (E(Fp2))';
                            1:Label4.Caption:='Public Key (E(Fp3))';
                            2:Label4.Caption:='Public Key (E(Fp4))';
                            end;
                            case Form3.ComboBox1.ItemIndex of
                            0:BLS1Generate;
                            1:BLS2Generate;
                            2:BLS3Generate;
                            end;
                            EnableInterface;
                           end;
end;

procedure TForm2.Button2Click(Sender: TObject);
begin
case ComboBox1.ItemIndex of
0:BLS1Sign;
1:BLS2Sign;
2:BLS3Sign;
end;
end;

procedure TForm2.Button3Click(Sender: TObject);
begin
case ComboBox1.ItemIndex of
0:BLS1Verify;
1:BLS2Verify;
2:BLS3Verify;
end;
end;

procedure TForm2.Button4Click(Sender: TObject);
var Output:Tbytes;
    btable1,btable2:Tbytes;
    Data:Byte;
    St:TStringList;
begin
Form4.CheckBox1.Show;
Form4.showmodal;
if (Form4.Modalresult=1)and(SaveDialog1.Execute) then
    begin
    if Form4.CheckBox1.Checked then Data :=1 shl 5 else Data:=0;
    Data:=Data Or (ComboBox1.ItemIndex or (ComboBox2.ItemIndex shl 2)) or(Byte(Form4.RadioButton2.Checked) shl 6);
    if Form4.CheckBox1.Checked then btable1:=HexToArray(Edit1.Text)
    else Setlength(btable1,0);
    btable2:=HexToArray(Edit3.Text);
    if Form4.RadioButton2.Checked then begin   // Binary Output format
                                          Setlength(Output,Length(btable1)+length(btable2)+1);
                                          if length(btable1)<>0 then Move(btable1[0],Output[1],length(btable1));
                                          Move(btable2[0],Output[Length(btable1)+1],Length(Btable2));
                                          Output[0]:=Data;
                                          TFile.WriteAllBytes(SaveDialog1.FileName,Output);
                                          end
    else begin   // Base-64 Output Format
         St:=TStringList.Create;
         St.Add('Algo :'+inttostr(Data));
         case ComboBox1.ItemIndex of
         0:St.Add('BLS1 Scheme (Bonhe-Lynn-Shacham[2001])');
         1:St.Add('BLS2 Scheme (Boneh-Boyen with RO[2004])');
         2:St.Add('BLS3 Scheme (Boneh-Boyen without RO[2004])');
         end;
         case ComboBox2.ItemIndex of
         0:St.Add('Security Level :128bit (BN Curve)');
         1:St.Add('Security Level :192bit (KSS18 Curve)');
         2:St.Add('Security Level :256bit (BLS24 Curve)');
         end;
         if length(btable1)<>0 then begin
                                    St.Add('Secret Key :');
                                    St.Add(EncodeBase64(btable1,length(btable1)));
                                    end;
         St.Add('Public Key :');
         St.Add(EncodeBase64(btable2,length(btable2)));
         St.SaveToFile(SaveDialog1.FileName);
         end;
    end;
end;

procedure TForm2.Button5Click(Sender: TObject);
var Input:Tbytes;
    btable1,btable2,tmp:Tbytes;
    Data:Byte;
    SecretSize,sbit,pbit,i:Integer;
    St:TStringList;
    Format,stmp:string;
begin
if OpenDialog1.Execute then begin
                            Input:=TFile.ReadAllBytes(OpenDialog1.Filename);
                            Data:=Input[0];
                            if Data =Ord('A') then begin   // Base-64 File Format
                                                   format:='b64';
                                                   St:=TStringList.Create;
                                                   St.LoadFromFile(OpenDialog1.FileName);
                                                   Try
                                                    Data:=Strtoint(Copy(St.Strings[0],Pos(':',St.Strings[0])+1,length(St.Strings[0])));
                                                    Except
                                                    MessageDlg('Format de fichier Invalide !',mtwarning,[mbok],0);
                                                    exit;
                                                    end;
                                                   end
                            else format:='bin';
                            ComboBox1.ItemIndex:=Data and 3;
                            ComboBox2.ItemIndex:=(Data shr 2) and 3;
                            case ComboBox1.ItemIndex of
                              0:RichEdit1.Text:='Boneh, Dan, Ben Lynn, and Hovav Shacham. "Short signatures from the Weil pairing." In International Conference on the Theory and Application of Cryptology and Information Security, pp. 514-532. '+'Springer, Berlin, Heidelberg, 2001. https://link.springer.com/article/10.1007/s00145-004-0314-9';
                              1:RichEdit1.Text:='Boneh, Dan, and Xavier Boyen. "Short signatures without random oracles." In International Conference on the Theory and Applications of Cryptographic Techniques, pp. 56-73. Springer, Berlin,'+' Heidelberg, 2004. https://link.springer.com/chapter/10.1007/978-3-540-24676-3_4';
                              2:RichEdit1.Text:='Boneh, Dan, and Xavier Boyen. "Short signatures without random oracles." In International Conference on the Theory and Applications of Cryptographic Techniques, pp. 56-73. Springer, Berlin,'+' Heidelberg, 2004. https://link.springer.com/chapter/10.1007/978-3-540-24676-3_4';
                            end;
                            Data:=Data shr 4;
                            if Data<>0 then begin    /// Secret Key Included
                                            case ComboBox2.ItemIndex of
                                             0:SecretSize:=BNCurvePairing1.CurveParams.R.BitLength div 8;
                                             1:SecretSize:=KSS18CurvePairing1.CurveParams.R.BitLength div 8;
                                             2:SecretSize:=BLS24CurvePairing1.CurveParams.R.BitLength div 8;
                                            end;
                                            While(SecretSize mod 8<>0) do SecretSize:=SecretSize+1;
                                            if ComboBox1.ItemIndex=2 then SecretSize:=2*SecretSize;
                                            Setlength(btable1,SecretSize);
                                            if format='bin' then Move(Input[1],btable1[0],Length(btable1))
                                            else begin
                                                 i:=4;
                                                 stmp:='';
                                                 While(Copy(St.Strings[i],1,6)<>'Public') do begin
                                                                                             stmp:=stmp+St.Strings[i];
                                                                                             i:=i+1;
                                                                                             end;
                                                 btable1:=DecodeBase64(stmp);
                                                 end;
                                            Edit1.Text:=ArrayToHex(btable1);
                                            end
                            else begin
                                 Edit1.Text:='';
                                 Setlength(btable1,0);
                                 i:=3;
                                 end;
                            Setlength(btable2,Length(Input)-Length(btable1)-1);
                            if format='bin' then Move(Input[length(btable1)+1],btable2[0],length(btable2))
                            else begin
                                 stmp:='';
                                 i:=i+1;
                                 while (i<St.Count) do begin
                                                       Stmp:=stmp+St.Strings[i];
                                                       i:=i+1;
                                                       end;
                                 btable2:=DecodeBase64(stmp);
                                 end;
                            Edit3.Text:=ArrayToHex(btable2);
                            case ComboBox2.ItemIndex of
                              0:Label4.Caption:='Public Key (E(Fp2))';
                              1:Label4.Caption:='Public Key (E(Fp3))';
                              2:Label4.Caption:='Public Key (E(Fp4))';
                            end;
                            if ComboBox1.ItemIndex<2 then SecretKey:=ByteArrayToLint(btable1)
                            else begin   // Splite Secret Key info to x,y if BLS3 is Userd
                                 Setlength(Tmp,Length(btable1) div 2);
                                 Move(btable1[0],Tmp[0],length(Tmp));
                                 SecretKey_BLSWithoutRO.x:=ByteArrayToLint(tmp);
                                 Move(btable1[Length(tmp)],Tmp[0],length(Tmp));
                                 SecretKey_BLSWithoutRO.y:=ByteArrayToLint(tmp);
                                 end;
                            case ComboBox2.ItemIndex of
                              0:begin   /// 128bit level : BN Curve
                                sbit:=BNCurvePairing1.CurveParams.R.BitLength;
                                pbit:=BNCurvePairing1.CurveParams.P.BitLength;
                                case ComboBox1.ItemIndex of
                                 0,1:begin //  BLS 1 / BLS 2 Signature Algo
                                     PublicKey_128BN.SetCurveParams(BNCurvePairing1.CurveParams);
                                     PublicKey_128BN.DeCompressFromArray(btable2);
                                     Edit2.Text:=PublicKey_128BN.ToHexString;
                                     Label6.Caption:=inttostr(sbit)+' bit';
                                     Label7.Caption:=inttostr(4*pbit)+' bit';
                                     Label8.Caption:=inttostr(2*pbit)+' bit';
                                     if ComboBox1.ItemIndex=1 then
                                      GeneratorsPairingsForBLS2_BN:=BNCurvePairing1.Paire(DefaultG1_128BN,DefaultG2_128BN);
                                     end;
                              2:begin  // BLS 3 Signature Algo
                                  PublicKey_BLSWithoutRO.Seclevel:=0;
                                  PublicKey_BLSWithoutRO.u_Fp2.SetCurveParams(BNCurvePairing1.CurveParams);
                                  PublicKey_BLSWithoutRO.v_Fp2.SetCurveParams(BNCurvePairing1.CurveParams);
                                  Setlength(Tmp,Length(btable2) div 2);
                                  Move(btable2[0],Tmp[0],length(Tmp));
                                  PublicKey_BLSWithoutRO.u_Fp2.DeCompressFromArray(tmp);
                                  Move(btable2[Length(tmp)],Tmp[0],length(Tmp));
                                  PublicKey_BLSWithoutRO.v_Fp2.DeCompressFromArray(tmp);
                                  Edit2.Text:='('+PublicKey_BLSWithoutRO.u_Fp2.ToHexString+' , '+PublicKey_BLSWithoutRO.v_Fp2.ToHexString+')';
                                  Label6.Caption:=inttostr(2*sbit)+' bit';
                                  Label7.Caption:=inttostr(8*pbit)+' bit';
                                  Label8.Caption:=inttostr(4*pbit)+' bit';
                                end;
                               end;
                               end;
                             1:begin    // 192bit level :KSS18 Curve
                               sbit:=KSS18CurvePairing1.CurveParams.R.BitLength;
                               pbit:=KSS18CurvePairing1.CurveParams.P.BitLength;
                               case ComboBox1.ItemIndex of
                                0,1:begin //  BLS 1 / BLS 2 Signature Algo
                                    PublicKey_192KSS.SetCurveParams(KSS18CurvePairing1.CurveParams);
                                    PublicKey_192KSS.DeCompressFromArray(btable2);
                                    Edit2.Text:=PublicKey_192KSS.ToHexString;
                                    Label6.Caption:=inttostr(sbit)+' bit';
                                    Label7.Caption:=inttostr(4*pbit)+' bit';
                                    Label8.Caption:=inttostr(2*pbit)+' bit';
                                    if ComboBox1.ItemIndex=1 then
                                     GeneratorsPairingsForBLS2_KSS:=KSS18CurvePairing1.Paire(DefaultG1_192KSS,DefaultG2_192KSS);
                                    end;
                               2:begin  // BLS 3 Signature Algo
                                   PublicKey_BLSWithoutRO.Seclevel:=1;
                                   PublicKey_BLSWithoutRO.u_Fp3.SetCurveParams(KSS18CurvePairing1.CurveParams);
                                   PublicKey_BLSWithoutRO.v_Fp3.SetCurveParams(KSS18CurvePairing1.CurveParams);
                                   Setlength(Tmp,Length(btable2) div 2);
                                   Move(btable2[0],Tmp[0],length(Tmp));
                                   PublicKey_BLSWithoutRO.u_Fp3.DeCompressFromArray(tmp);
                                   Move(btable2[Length(tmp)],Tmp[0],length(Tmp));
                                   PublicKey_BLSWithoutRO.v_Fp3.DeCompressFromArray(tmp);
                                   Edit2.Text:='('+PublicKey_BLSWithoutRO.u_Fp3.ToHexString+' , '+PublicKey_BLSWithoutRO.v_Fp3.ToHexString+')';
                                   Label6.Caption:=inttostr(2*sbit)+' bit';
                                   Label7.Caption:=inttostr(8*pbit)+' bit';
                                   Label8.Caption:=inttostr(4*pbit)+' bit';
                                 end;
                                end;
                                end;
                             2:begin   // 256bit level : BLS24 Curve
                               sbit:=BLS24CurvePairing1.CurveParams.R.BitLength;
                               pbit:=BLS24CurvePairing1.CurveParams.P.BitLength;
                               case ComboBox1.ItemIndex of
                                0,1:begin //  BLS 1 /BLS 2 Signature Algo
                                    PublicKey_256NBLS.SetCurveParams(BLS24CurvePairing1.CurveParams);
                                    PublicKey_256NBLS.DeCompressFromArray(btable2);
                                    Edit2.Text:=PublicKey_256NBLS.ToHexString;
                                    Label6.Caption:=inttostr(2*sbit)+' bit';
                                    Label7.Caption:=inttostr(8*pbit)+' bit';
                                    Label8.Caption:=inttostr(4*pbit)+' bit';
                                    if ComboBox1.ItemIndex=1 then
                                     GeneratorsPairingsForBLS2_BLS24:=BLS24CurvePairing1.Paire(DefaultG1_256BLS,DefaultG2_256BLS);
                                    end;
                                2:begin  // BLS 3 Signature Algo
                                    PublicKey_BLSWithoutRO.Seclevel:=2;
                                    PublicKey_BLSWithoutRO.u_Fp4.SetCurveParams(BLS24CurvePairing1.CurveParams);
                                    PublicKey_BLSWithoutRO.v_Fp4.SetCurveParams(BLS24CurvePairing1.CurveParams);
                                    Setlength(Tmp,Length(btable2) div 2);
                                    Move(btable2[0],Tmp[0],length(Tmp));
                                    PublicKey_BLSWithoutRO.u_Fp4.DeCompressFromArray(tmp);
                                    Move(btable2[Length(tmp)],Tmp[0],length(Tmp));
                                    PublicKey_BLSWithoutRO.v_Fp4.DeCompressFromArray(tmp);
                                    Edit2.Text:='('+PublicKey_BLSWithoutRO.u_Fp4.ToHexString+' , '+PublicKey_BLSWithoutRO.v_Fp4.ToHexString+')';
                                    Label6.Caption:=inttostr(2*sbit)+' bit';
                                    Label7.Caption:=inttostr(8*pbit)+' bit';
                                    Label8.Caption:=inttostr(4*pbit)+' bit';
                                  end;
                                 end;
                                 end;
                                end;
                               EnableInterface;
                               end;
end;

procedure TForm2.Button6Click(Sender: TObject);
begin
case ComboBox2.ItemIndex of
0:form1.TreeView11.Curve:=BNCurvePairing1;
1:form1.TreeView11.Curve:=KSS18CurvePairing1;
2:form1.TreeView11.Curve:=BLS24CurvePairing1;
end;
Form1.ShowModal;
end;

procedure TForm2.EnableInterface;
begin
Button4.Enabled:=true;
Button2.Enabled:=(Edit1.Text<>'');
Button3.Enabled:=(Edit2.Text<>'');
Label12.Enabled:=true;
Label13.Enabled:=true;
Label3.Enabled:=true;
Label4.Enabled:=true;
Label5.Enabled:=true;
Label6.Enabled:=true;
Label7.Enabled:=true;
Label8.Enabled:=true;
Label9.Enabled:=true;
Label10.Enabled:=true;
Label11.Enabled:=true;
Memo1.Enabled:=true;
Button7.Enabled:=true;
Edit4.Enabled:=true;
SpeedButton1.Enabled:=true;
SpeedButton2.Enabled:=true;
end;

procedure TForm2.FormCreate(Sender: TObject);
var
  mask : Integer;
begin
  mask := SendMessage(RichEdit1.Handle, EM_GETEVENTMASK, 0, 0);
  SendMessage(RichEdit1.Handle, EM_SETEVENTMASK, 0, mask or ENM_LINK);
  SendMessage(RichEdit1.Handle, EM_AUTOURLDETECT, Integer(True), 0);
//
DefaultG1_128BN.SetCurveParams(BNCurvePairing1.CurveParams);
DefaultG1_128BN.SetToDefaultGenerator;
DefaultG2_128BN.SetCurveParams(BNCurvePairing1.CurveParams);
DefaultG2_128BN.SetToDefaultGenerator;

DefaultG1_192KSS.SetCurveParams(KSS18CurvePairing1.CurveParams);
DefaultG1_192KSS.SetToDefaultGenerator;
DefaultG2_192KSS.SetCurveParams(KSS18CurvePairing1.CurveParams);
DefaultG2_192KSS.SetToDefaultGenerator;

DefaultG1_256BLS.SetCurveParams(BLS24CurvePairing1.CurveParams);
DefaultG1_256BLS.SetToDefaultGenerator;
DefaultG2_256BLS.SetCurveParams(BLS24CurvePairing1.CurveParams);
DefaultG2_256BLS.SetToDefaultGenerator;

//
end;
procedure TForm2.SpeedButton1Click(Sender: TObject);
Var st:TStringList;
    stmp:string;
    btable:Tbytes;
begin
if OpenDialog1.Execute then begin
                            stmp:=TFile.ReadAllText(OpenDialog1.FileName);
                            if copy(stmp,1,3)='BLS' then begin
                                                         St:=TStringList.Create;
                                                         St.LoadFromFile(SaveDialog1.FileName);
                                                         Edit4.Text:=ArrayToHex(DecodeBase64(St.Strings[St.Count-1]));
                                                         end
                            else begin
                                 btable:=TFile.ReadAllBytes(OpenDialog1.FileName);
                                 Edit4.Text:=ArrayToHex(btable);
                                 end;

                            end;
end;

procedure TForm2.SpeedButton2Click(Sender: TObject);
var St: TStringList;
    btable:Tbytes;
begin
Form4.CheckBox1.Hide;
if (Form4.ShowModal=1)and(SaveDialog1.Execute()) then begin
                          if Form4.RadioButton2.Checked then TFile.WriteAllBytes(SaveDialog1.FileName,HexToArray(Edit4.Text))
                          else begin
                               st:=TStringList.Create;
                               case ComboBox1.ItemIndex of
                               0:St.Add('BLS1 Scheme (Bonhe-Lynn-Shacham[2001])');
                               1:St.Add('BLS2 Scheme (Boneh-Boyen with RO[2004])');
                               2:St.Add('BLS3 Scheme (Boneh-Boyen without RO[2004])');
                               end;
                               case ComboBox2.ItemIndex of
                               0:St.Add('Security Level :128bit (BN Curve)');
                               1:St.Add('Security Level :192bit (KSS18 Curve)');
                               2:St.Add('Security Level :256bit (BLS24 Curve)');
                               end;
                               St.Add('BLS Signature:');
                               btable:=HexToArray(Edit4.Text);
                               St.Add(EncodeBase64(btable,length(btable)));
                               St.SaveToFile(SaveDialog1.FileName);
                               end;
                          end;
end;

procedure TForm2.Split(const Delimiter: Char; Input: string; const Strings: TStrings) ;
begin
   Assert(Assigned(Strings)) ;
   Strings.Clear;
   Strings.Delimiter := Delimiter;
   Strings.DelimitedText := Input;
end;
procedure Tform2.WndProc(var Msgs: TMessage);
var
  p: TENLink;
  strURL, link, path : string;
  A : TStringList;
begin
  if (Msgs.Msg = WM_NOTIFY) then
  begin
    if (PNMHDR(Msgs.LParam).code = EN_LINK) then
    begin
      p := TENLink(Pointer(TWMNotify(Msgs).NMHdr)^);
      if (p.msg = WM_LBUTTONDOWN) then
      begin
        SendMessage(RichEdit1.Handle, EM_EXSETSEL, 0, LongInt(@(p.chrg)));
        link := RichEdit1.SelText;
        A := TStringList.Create;
        Split(':', link, A) ;
        path := GetEnvironmentVariable('TEMP');
        strURL := 'file:\\' + path + '\' + A[1] + '.png';
        A.Free;
        ShowMessage(strURL);
        ShellExecute(Handle, 'open', PChar(strURL), 0, 0, SW_SHOWNORMAL);
      end
    end
  end;

  inherited;
end;

//*******************************************************************************************************************//
// First Scheme of BLS Signature
//  Boneh, Dan, Ben Lynn, and Hovav Shacham. "Short signatures from the Weil pairing."
//*******************************************************************************************************************//
          // ***************************** Key Pair Generation ************************///
Procedure TForm2.BLS1Generate;
Var btable:TBytes;
    sbit,pbit:Integer;
begin
case ComboBox2.ItemIndex of
0:begin
  SecretKey.SetToRandom(BNCurvePairing1.CurveParams.R);
  PublicKey_128BN:=SecretKey*DefaultG2_128BN;
  Edit2.Text:=PublicKey_128BN.ToHexString;
  btable:=PublicKey_128BN.CompressToArray;
  sbit:=BNCurvePairing1.CurveParams.R.BitLength;
  pbit:=BNCurvePairing1.CurveParams.P.BitLength;
  end;
1:begin
  SecretKey.SetToRandom(KSS18CurvePairing1.CurveParams.R);
  PublicKey_192KSS:=SecretKey*DefaultG2_192KSS;
  Edit2.Text:=PublicKey_192KSS.ToHexString;
  btable:=PublicKey_192KSS.CompressToArray;
  sbit:=KSS18CurvePairing1.CurveParams.R.BitLength;
  pbit:=KSS18CurvePairing1.CurveParams.P.BitLength;
  end;
2:begin
  SecretKey.SetToRandom(BLS24CurvePairing1.CurveParams.R);
  PublicKey_256NBLS:=SecretKey*DefaultG2_256BLS;
  Edit2.Text:=PublicKey_256NBLS.ToHexString;
  btable:=PublicKey_256NBLS.CompressToArray;
  sbit:=BLS24CurvePairing1.CurveParams.R.BitLength;
  pbit:=BLS24CurvePairing1.CurveParams.P.BitLength;
  end;
end;
Edit1.Text:=SecretKey.ToHexString;
Edit3.Text:= ArrayToHex(btable);
Label6.Caption:=inttostr(sbit)+' bit';
Label7.Caption:=inttostr(4*pbit)+' bit';
Label8.Caption:=inttostr(2*pbit)+' bit';
end;
        // ***************************** BLS1 Signature ************************///
Procedure TForm2.BLS1Sign;
var MessageHash,FileSignature:TBytes;
    SigAsPoint,tmp:G1BN;
begin
tim:=now;
MessageHash:=SHA256StringHash(Memo1.Text);
case ComboBox2.ItemIndex of
0:tmp.SetCurveParams(BNCurvePairing1.CurveParams);
1:tmp.SetCurveParams(KSS18CurvePairing1.CurveParams);
2:tmp.SetCurveParams(BLS24CurvePairing1.CurveParams);
end;
tmp.SetAsTorsionFromHash(MessageHash);
//SigAsPoint:=SecretKey*SigAsPoint;
_Mul_NAF_Fp_FpPoint (SecretKey,tmp,SigAsPoint);
FileSignature:=SigAsPoint.CompressToArray;
Edit4.Text:=ArrayToHex(FileSignature);
tim:=now-tim;
Label11.Caption:=inttostr((length(FileSignature)-1)*8)+' bit';
MessageDlg('Signature Générée avec succés ! (en '+TimetoMessage(Tim)+')',mtConfirmation,[mbok],0);
end;

        // ***************************** BLS1 Verification ************************///
Procedure TForm2.BLS1Verify;
var FileSignature,Filehash:TBytes;
    SigAsPoint:G1BN;
    hashAspoint:G1BN;
    Equals:boolean;
    e1_12,e2_12:GTBN;
    e1_18,e2_18:GTKSS18;
    e1_24,e2_24:GtBls24;
begin
Tim:=now;
case ComboBox2.ItemIndex of
  0:begin
    SigAsPoint.SetCurveParams(BNCurvePairing1.CurveParams);
    hashAspoint.SetCurveParams(BNCurvePairing1.CurveParams);
    end;
  1:begin
    SigAsPoint.SetCurveParams(KSS18CurvePairing1.CurveParams);
    hashAspoint.SetCurveParams(KSS18CurvePairing1.CurveParams);
    end;
  2:begin
    SigAsPoint.SetCurveParams(BLS24CurvePairing1.CurveParams);
    hashAspoint.SetCurveParams(BLS24CurvePairing1.CurveParams);
    end;
  end;
  FileHash:=SHA256StringHash(Memo1.Text);
  hashAspoint.SetAsTorsionFromHash(FileHash);
  FileSignature:=HexToArray(Edit4.Text);
  SigAsPoint.DeCompressFromArray(FileSignature);
case ComboBox2.ItemIndex of
  0:begin
    e1_12:=BNCurvePairing1.Paire(SigAsPoint,DefaultG2_128BN);
    e2_12:=BNCurvePairing1.Paire(hashAspoint,PublicKey_128BN);
    Equals:=e1_12=e2_12;
    end;
  1:begin
    e1_18:=KSS18CurvePairing1.Paire(SigAsPoint,DefaultG2_192KSS);
    e2_18:=KSS18CurvePairing1.Paire(hashAspoint,PublicKey_192KSS);
    Equals:=e1_18=e2_18;
    end;
  2:begin
    e1_24:=BLS24CurvePairing1.Paire(SigAsPoint,DefaultG2_256BLS);
    e2_24:=BLS24CurvePairing1.Paire(hashAspoint,PublicKey_256NBLS);
    Equals:=e1_24=e2_24;
    end;
  end;
Tim:=now-Tim;
if  Equals then MessageDlg('Signature Valide :'+#13+#10+ArrayToHex(FileSignature)+#13+#10+'Vérification faite en : '+TimetoMessage(Tim)+'.',mtInformation,[mbok],0)
else Messagedlg('Signature Invalide :la signature ne correspond pas à la clé de vérification, ou bien au message signé',mterror,[mbok],0);
end;

//*******************************************************************************************************************//
// Second Scheme of BLS Signature  " With RO "
// Boneh, Dan, and Xavier Boyen. "Short signatures without random oracles."
//*******************************************************************************************************************//
          // ***************************** Key Pair Generation ************************///
Procedure TForm2.BLS2Generate;
Var btable:TBytes;
    sbit,pbit:Integer;
begin
case ComboBox2.ItemIndex of
0:begin
  SecretKey.SetToRandom(BNCurvePairing1.CurveParams.R);
  PublicKey_128BN:=SecretKey*DefaultG2_128BN;
  Edit2.Text:=PublicKey_128BN.ToHexString;
  btable:=PublicKey_128BN.CompressToArray;
  sbit:=BNCurvePairing1.CurveParams.R.BitLength;
  pbit:=BNCurvePairing1.CurveParams.P.BitLength;
  GeneratorsPairingsForBLS2_BN:=BNCurvePairing1.Paire(DefaultG1_128BN,DefaultG2_128BN);
  end;
1:begin
  SecretKey.SetToRandom(KSS18CurvePairing1.CurveParams.R);
  PublicKey_192KSS:=SecretKey*DefaultG2_192KSS;
  Edit2.Text:=PublicKey_192KSS.ToHexString;
  btable:=PublicKey_192KSS.CompressToArray;
  sbit:=KSS18CurvePairing1.CurveParams.R.BitLength;
  pbit:=KSS18CurvePairing1.CurveParams.P.BitLength;
  GeneratorsPairingsForBLS2_KSS:=KSS18CurvePairing1.Paire(DefaultG1_192KSS,DefaultG2_192KSS);
  end;
2:begin
  SecretKey.SetToRandom(BLS24CurvePairing1.CurveParams.R);
  PublicKey_256NBLS:=SecretKey*DefaultG2_256BLS;
  Edit2.Text:=PublicKey_256NBLS.ToHexString;
  btable:=PublicKey_256NBLS.CompressToArray;
  sbit:=BLS24CurvePairing1.CurveParams.R.BitLength;
  pbit:=BLS24CurvePairing1.CurveParams.P.BitLength;
  GeneratorsPairingsForBLS2_BLS24:=BLS24CurvePairing1.Paire(DefaultG1_256BLS,DefaultG2_256BLS);
  end;
end;
Edit1.Text:=SecretKey.ToHexString;
Edit3.Text:= ArrayToHex(btable);
Label6.Caption:=inttostr(sbit)+' bit';
Label7.Caption:=inttostr(4*pbit)+' bit';
Label8.Caption:=inttostr(2*pbit)+' bit';
end;
               // ***************************** BLS2 Signature ************************///
Procedure TForm2.BLS2Sign;
var MessageHash,MessageSignature:TBytes;
    SigAsPoint:G1BN;
    m:LargeInt;
begin
Tim:=now;
MessageHash:=SHA256StringHash(Memo1.Text);
case ComboBox2.ItemIndex of
0:begin
  SigAsPoint.SetCurveParams(BNCurvePairing1.CurveParams);
  m:=(ByteArrayToLint(MessageHash)mod SigAsPoint.CurveParams.R);
  m:=(m+SecretKey) mod SigAsPoint.CurveParams.R;
  m:=m.InversModulo(SigAsPoint.CurveParams.R);
  _Mul_NAF_Fp_FpPoint(m,DefaultG1_128BN,SigAsPoint);
  end;
1:begin
  SigAsPoint.SetCurveParams(KSS18CurvePairing1.CurveParams);
  m:=(ByteArrayToLint(MessageHash) mod SigAsPoint.CurveParams.R)+SecretKey;
  m:=m.InversModulo(SigAsPoint.CurveParams.R);
  _Mul_NAF_Fp_FpPoint(m,DefaultG1_192KSS,SigAsPoint);
  end;
2:begin
  SigAsPoint.SetCurveParams(BLS24CurvePairing1.CurveParams);
  m:=(ByteArrayToLint(MessageHash) mod SigAsPoint.CurveParams.R)+SecretKey;
  m:=m.InversModulo(SigAsPoint.CurveParams.R);
  _Mul_NAF_Fp_FpPoint(m,DefaultG1_256BLS,SigAsPoint);
  end;
end;
MessageSignature:=SigAsPoint.CompressToArray;
Tim:=now-Tim;
Edit4.Text:=ArrayToHex(MessageSignature);
Label11.Caption:=inttostr((length(MessageSignature)-1)*8)+' bit';
MessageDlg('Signature Générée avec succés ! (en '+TimetoMessage(Tim)+')',mtConfirmation,[mbok],0);
end;

        // ***************************** BLS2 Verification ************************///
Procedure TForm2.BLS2Verify;
var FileSignature,Filehash:TBytes;
    SigAsPoint:G1BN;
    hashAsInt:LargeInt;
    Equals:boolean;
    e1_12,e2_12:GTBN;
    e1_18,e2_18:GTKSS18;
    e1_24,e2_24:GtBls24;
    TmpG2_BN:G2BN;
    TmpG2_KSS:G2KSS18;
    TmpG2_BLS:G2BLS24;
begin
Tim:=now;
FileHash:=SHA256StringHash(Memo1.Text);
case ComboBox2.ItemIndex of
  0:begin
    SigAsPoint.SetCurveParams(BNCurvePairing1.CurveParams);
    hashAsInt:=ByteArrayToLint(Filehash) mod BNCurvePairing1.CurveParams.R;
    _Mul_NAF_Fp2_FpPoint(hashAsInt,DefaultG2_128BN,TmpG2_BN,csProjective);
    TmpG2_BN:=TmpG2_BN+PublicKey_128BN;
    end;
  1:begin
    SigAsPoint.SetCurveParams(KSS18CurvePairing1.CurveParams);
    hashAsInt:=ByteArrayToLint(Filehash) mod KSS18CurvePairing1.CurveParams.R;
     _Mul_NAF_Fp3_FpPoint(hashAsInt,DefaultG2_192KSS,TmpG2_KSS,csProjective);
    TmpG2_KSS:=TmpG2_KSS+PublicKey_192KSS;
    end;
  2:begin
    SigAsPoint.SetCurveParams(BLS24CurvePairing1.CurveParams);
    hashAsInt:=ByteArrayToLint(Filehash) mod BLS24CurvePairing1.CurveParams.R;
    _Mul_NAF_Fp_Fp4Point(hashAsInt,DefaultG2_256BLS,TmpG2_BLS,csProjective);
    TmpG2_BLS:=TmpG2_BLS+PublicKey_256NBLS;
    end;
  end;
FileSignature:=HexToArray(Edit4.Text);
SigAsPoint.DeCompressFromArray(FileSignature);
case ComboBox2.ItemIndex of
  0:begin
    e1_12:=BNCurvePairing1.Paire(SigAsPoint,TmpG2_BN);
    Equals:=e1_12=GeneratorsPairingsForBLS2_BN;
    end;
  1:begin
    e1_18:=KSS18CurvePairing1.Paire(SigAsPoint,TmpG2_KSS);
    Equals:=e1_18=GeneratorsPairingsForBLS2_KSS;
    end;
  2:begin
    e1_24:=BLS24CurvePairing1.Paire(SigAsPoint,TmpG2_BLS);
    Equals:=e1_24=GeneratorsPairingsForBLS2_BLS24;
    end;
  end;
Tim:=now-tim;
if  Equals then MessageDlg('Signature Valide :'+#13+#10+ArrayToHex(FileSignature)+#13+#10+'Vérification faite en : '+TimetoMessage(Tim)+'.',mtInformation,[mbok],0)
else Messagedlg('Signature Invalide :la signature ne correspond pas à la clé de vérification, ou bien au message signé',mterror,[mbok],0);
end;

//*******************************************************************************************************************//
// Third Scheme of BLS Signature  "Without RO "
// Boneh, Dan, and Xavier Boyen. "Short signatures without random oracles."
//*******************************************************************************************************************//
          // ***************************** Key Pair Generation ************************///
Procedure TForm2.BLS3Generate;
Var btable1,btable2,merged:TBytes;
    sbit,pbit:Integer;
begin
case ComboBox2.ItemIndex of
0:begin
  SecretKey_BLSWithoutRO.x.SetToRandom(BNCurvePairing1.CurveParams.R);
  SecretKey_BLSWithoutRO.y.SetToRandom(BNCurvePairing1.CurveParams.R);
  PublicKey_BLSWithoutRO.Seclevel:=0;
  PublicKey_BLSWithoutRO.u_Fp2:=SecretKey_BLSWithoutRO.x*DefaultG2_128BN;
  PublicKey_BLSWithoutRO.v_Fp2:=SecretKey_BLSWithoutRO.y*DefaultG2_128BN;
  Edit2.Text:='('+PublicKey_BLSWithoutRO.u_Fp2.ToHexString+' , '+PublicKey_BLSWithoutRO.v_Fp2.ToHexString+')';
  btable1:=PublicKey_BLSWithoutRO.u_Fp2.CompressToArray;
  btable2:=PublicKey_BLSWithoutRO.v_Fp2.CompressToArray;
  Setlength(merged,length(btable1)+length(btable2));
  Move(btable1[0],merged[0],length(btable1));
  Move(btable2[0],merged[length(btable1)],length(btable2));
  sbit:=BNCurvePairing1.CurveParams.R.BitLength;
  pbit:=BNCurvePairing1.CurveParams.P.BitLength;
  GeneratorsPairingsForBLS2_BN:=BNCurvePairing1.Paire(DefaultG1_128BN,DefaultG2_128BN);
  end;
1:begin
  SecretKey_BLSWithoutRO.x.SetToRandom(KSS18CurvePairing1.CurveParams.R);
  SecretKey_BLSWithoutRO.y.SetToRandom(KSS18CurvePairing1.CurveParams.R);
  PublicKey_BLSWithoutRO.Seclevel:=1;
  PublicKey_BLSWithoutRO.u_Fp3:=SecretKey_BLSWithoutRO.x*DefaultG2_192KSS;
  PublicKey_BLSWithoutRO.v_Fp3:=SecretKey_BLSWithoutRO.y*DefaultG2_192KSS;
  Edit2.Text:='('+PublicKey_BLSWithoutRO.u_Fp3.ToHexString+' , '+PublicKey_BLSWithoutRO.v_Fp3.ToHexString+')';
  btable1:=PublicKey_BLSWithoutRO.u_Fp3.CompressToArray;
  btable2:=PublicKey_BLSWithoutRO.v_Fp3.CompressToArray;
  Setlength(merged,length(btable1)+length(btable2));
  Move(btable1[0],merged[0],length(btable1));
  Move(btable2[0],merged[length(btable1)],length(btable2));
  sbit:=KSS18CurvePairing1.CurveParams.R.BitLength;
  pbit:=KSS18CurvePairing1.CurveParams.P.BitLength;
  GeneratorsPairingsForBLS2_KSS:=KSS18CurvePairing1.Paire(DefaultG1_192KSS,DefaultG2_192KSS);
  end;
2:begin
  SecretKey_BLSWithoutRO.x.SetToRandom(BLS24CurvePairing1.CurveParams.R);
  SecretKey_BLSWithoutRO.y.SetToRandom(BLS24CurvePairing1.CurveParams.R);
  PublicKey_BLSWithoutRO.Seclevel:=2;
  PublicKey_BLSWithoutRO.u_Fp4:=SecretKey_BLSWithoutRO.x*DefaultG2_256BLS;
  PublicKey_BLSWithoutRO.v_Fp4:=SecretKey_BLSWithoutRO.y*DefaultG2_256BLS;
  Edit2.Text:='('+PublicKey_BLSWithoutRO.u_Fp4.ToHexString+' , '+PublicKey_BLSWithoutRO.v_Fp4.ToHexString+')';
  btable1:=PublicKey_BLSWithoutRO.u_Fp4.CompressToArray;
  btable2:=PublicKey_BLSWithoutRO.v_Fp4.CompressToArray;
  Setlength(merged,length(btable1)+length(btable2));
  Move(btable1[0],merged[0],length(btable1));
  Move(btable2[0],merged[length(btable1)],length(btable2));
  sbit:=BLS24CurvePairing1.CurveParams.R.BitLength;
  pbit:=BLS24CurvePairing1.CurveParams.P.BitLength;
  GeneratorsPairingsForBLS2_BLS24:=BLS24CurvePairing1.Paire(DefaultG1_256BLS,DefaultG2_256BLS);
  end;
end;
Edit3.Text:= ArrayToHex(merged);
btable1:=LIntToByteArray(SecretKey_BLSWithoutRO.x);
btable2:=LIntToByteArray(SecretKey_BLSWithoutRO.y);
Setlength(merged,length(btable1)+length(btable2));
Move(btable1[0],merged[0],length(btable1));
Move(btable2[0],merged[length(btable1)],length(btable2));
Edit1.Text:=ArrayToHex(merged);
Label6.Caption:=inttostr(2*sbit)+' bit';
Label7.Caption:=inttostr(8*pbit)+' bit';
Label8.Caption:=inttostr(4*pbit)+' bit';
end;

               // ***************************** BLS3 Signature ************************///
procedure TForm2.BLS3Sign;
var MessageHash,MessageSignature:TBytes;
    SigAsPoint:G1BN;
    m,r:LargeInt;
    btable1,btable2,merged:TBytes;
begin
tim:=now;
MessageHash:=SHA256StringHash(Memo1.Text);
case ComboBox2.ItemIndex of
0:begin
  SigAsPoint.SetCurveParams(BNCurvePairing1.CurveParams);
  r.SetToRandom(SigAsPoint.CurveParams.R);
  m:=(ByteArrayToLint(MessageHash)+r*SecretKey_BLSWithoutRO.y+SecretKey_BLSWithoutRO.x) mod BNCurvePairing1.CurveParams.R;
  m:=m.InversModulo(SigAsPoint.CurveParams.R);
  SigAsPoint:=m*DefaultG1_128BN;
  end;
1:begin
  SigAsPoint.SetCurveParams(KSS18CurvePairing1.CurveParams);
  r.SetToRandom(SigAsPoint.CurveParams.R);
  m:=(ByteArrayToLint(MessageHash)+r*SecretKey_BLSWithoutRO.y+SecretKey_BLSWithoutRO.x)mod KSS18CurvePairing1.CurveParams.R;
  m:=m.InversModulo(SigAsPoint.CurveParams.R);
  SigAsPoint:=m*DefaultG1_192KSS;
  end;
2:begin
  SigAsPoint.SetCurveParams(BLS24CurvePairing1.CurveParams);
  r.SetToRandom(SigAsPoint.CurveParams.R);
  m:=(ByteArrayToLint(MessageHash)+r*SecretKey_BLSWithoutRO.y+SecretKey_BLSWithoutRO.x) mod BLS24CurvePairing1.CurveParams.R;
  m:=m.InversModulo(SigAsPoint.CurveParams.R);
  SigAsPoint:=m*DefaultG1_256BLS;
  end;
end;
btable1:=SigAsPoint.CompressToArray;
btable2:=LIntToByteArray(r);
Setlength(MessageSignature,length(btable1)+length(btable2));
Move(btable2[0],MessageSignature[0],length(btable2));
Move(btable1[0],MessageSignature[length(btable2)],length(btable1));
tim:=now-tim;
Edit4.Text:=ArrayToHex(MessageSignature);
Label11.Caption:=inttostr((length(MessageSignature)-1)*8)+' bit';
MessageDlg('Signature Générée avec succés ! (en '+TimetoMessage(Tim)+')',mtConfirmation,[mbok],0);
end;
        // ***************************** BLS3 Verification ************************///
Procedure TForm2.BLS3Verify;
var FileSignature,Filehash:TBytes;
    btable1,btable2,merged:TBytes;
    SigAsPoint:G1BN;
    hashAsInt,R:LargeInt;
    Rsize:Word;
    Equals:boolean;
    e1_12,e2_12:GTBN;
    e1_18,e2_18:GTKSS18;
    e1_24,e2_24:GtBls24;
    TmpG2:BLSwithoutROPublicKey;
begin
Tim:=now;
FileHash:=SHA256StringHash(Memo1.Text);
case ComboBox2.ItemIndex of
  0:SigAsPoint.SetCurveParams(BNCurvePairing1.CurveParams);
  1:SigAsPoint.SetCurveParams(KSS18CurvePairing1.CurveParams);
  2:SigAsPoint.SetCurveParams(BLS24CurvePairing1.CurveParams);
  end;
hashAsInt:=ByteArrayToLint(Filehash) mod SigAsPoint.CurveParams.R;
Rsize:=SigAsPoint.CurveParams.R.BitLength;
while (Rsize mod 8<>0) do Rsize:=Rsize+1;
merged:=HexToArray(Edit4.Text);
Setlength(btable1,Rsize div 8);
Setlength(btable2,length(merged) -(Rsize div 8));
Move(merged[0],btable1[0],Rsize div 8);
Move(merged[Rsize div 8],btable2[0],Length(btable2));
SigAsPoint.DeCompressFromArray(btable2);
R:=ByteArrayToLint(btable1);
case ComboBox2.ItemIndex of
  0:begin
    _Mul_Fp_Fp2Point(R,PublicKey_BLSWithoutRO.v_Fp2,TmpG2.v_Fp2);
    _Mul_Fp_Fp2Point(hashAsInt,DefaultG2_128BN,TmpG2.u_Fp2);
    _Add_Fp2_Point(TmpG2.u_Fp2,TmpG2.v_Fp2,TmpG2.u_Fp2);
    _Add_Fp2_Point(TmpG2.u_Fp2,PublicKey_BLSWithoutRO.u_Fp2,TmpG2.u_Fp2);
    e1_12:=BNCurvePairing1.Paire(SigAsPoint,TmpG2.u_Fp2);
    Equals:=e1_12=GeneratorsPairingsForBLS2_BN;
    end;
  1:begin
    _Mul_Fp_Fp3Point(R,PublicKey_BLSWithoutRO.v_Fp3,TmpG2.v_Fp3);
    _Mul_Fp_Fp3Point(hashAsInt,DefaultG2_192KSS,TmpG2.u_Fp3);
    _Add_Fp3_Point(TmpG2.u_Fp3,TmpG2.v_Fp3,TmpG2.u_Fp3);
    _Add_Fp3_Point(TmpG2.u_Fp3,PublicKey_BLSWithoutRO.u_Fp3,TmpG2.u_Fp3);
    e1_18:=KSS18CurvePairing1.Paire(SigAsPoint,TmpG2.u_Fp3);
    Equals:=e1_18=GeneratorsPairingsForBLS2_KSS;
    end;
  2:begin
    _Mul_Fp_Fp4Point(R,PublicKey_BLSWithoutRO.v_Fp4,TmpG2.v_Fp4);
    _Mul_Fp_Fp4Point(hashAsInt,DefaultG2_256BLS,TmpG2.u_Fp4);
    _Add_Fp4_Point(TmpG2.u_Fp4,TmpG2.v_Fp4,TmpG2.u_Fp4);
    _Add_Fp4_Point(TmpG2.u_Fp4,PublicKey_BLSWithoutRO.u_Fp4,TmpG2.u_Fp4);
    e1_24:=BLS24CurvePairing1.Paire(SigAsPoint,TmpG2.u_Fp4);
    Equals:=e1_24=GeneratorsPairingsForBLS2_BLS24;
    end;
  end;
Tim:=now-tim;
if  Equals then MessageDlg('Signature Valide :'+#13+#10+ArrayToHex(merged)+#13+#10+'Vérification faite en : '+TimetoMessage(Tim)+'.',mtInformation,[mbok],0)
else Messagedlg('Signature Invalide :la signature ne correspond pas à la clé de vérification, ou bien au message signé',mterror,[mbok],0);
end;


end.
