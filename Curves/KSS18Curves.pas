unit KSS18Curves;

interface

uses  Vcl.ComCtrls, System.SysUtils,VLargeIntegers,LargeIntegers,Fp18Arithmetic,Fp9Arithmetic,Fp3Arithmetic,
      GeneralTypes,System.classes,VCL.dialogs,ECCFp,Eccfp3,TreeView1;


Type
   GTKSS18=Fp18Int;
   G2KSS18=Fp3Point;
   G1KSS18=FpPoint;
   LargeInt=Lint;
   KSS18CurvesPairingAlgos=(KSS18pOptAte);
   StandardKSS18Curves=(scKSS18at256,scKSS18at192_1,scKSS18at192_2,scKSS18at192_3,scKSS18at128);

   KSS18CurvesParamsDefinition=record
                            SecurityLevel:String;
                            u:String;  // the paramater of generation for the KSS18 curve
                            Beta:Integer; // non-square elements of the irredictible polynomial on FP2
                            Sigma:string; // non-square non-cube elements of the irredictible polynomial on FP4
                            Gamma:string; // non-square non-cube elements of the irredictible polynomial on FP8
                            A,B:String; // parametres of the curve   y^2=x^3+Ax+B
                            TwistMode:TTwistModel;
                            GenratorX:String;
                            TwistGeneratorSeed:Word;
                            end;
   ListOfKSS18Params=array of KSS18CurvesParamsDefinition;

 TKSS18CurvePairing=class (TCurve)
                      {Definition of a Curve with respect the the Weistrass model   Y^2=X^3+A*X+B (mod P)}
              private
              _A,_B:Lint;
              _u:Lint;
              _LoopMode:LoopPoweringMode;                     //  Loopmode :Negative/Binary Representation
              _FP18PoweringMode:FpPoweringMode;             // Final Exponentiation Powering Mode : Beuchat Approach, Karabina Approach, Classical Square and Multiply
              _PairingAlgorithm:KSS18CurvesPairingAlgos;         //  Pairing Algorithme OptAte,R-Ate.....
              Loop:PLIntArrayForm;
              params:StandardKSS18Curves;
              procedure SetA(Value:String);
              function getA:String;
              procedure SetB(Value:String);
              function getB:String;
              procedure Setu(Value:String);
              function getu:String;
              function getRtw:String;
              function getAtw:String;
              function getBtw:String;
              procedure SetSigma(Value:String);
              function getSigma:String;
              function getGamma:String;
              procedure SetGamma(Value:String);
              procedure SetBeta(Value:integer);
              function getBeta:integer;
              function getR:string;
              function getN:string;
              function getTr:string;
              function getLp:string;
              function getHtw:string;
              function getH:string;
              function getP:string;
              procedure setCoord(value:CoordinatesSystem);
              function getCoord:CoordinatesSystem;
              function GetTwm:TTwistModel;
              procedure RecomputeParametres;
              Procedure FinalPowerOptimalAteKSS18(ff:Fp18Int;u:LInt;PoweringMode:FpPoweringMode;var Result:Fp18Int); // Compute f^((p^18-1)/r)
              Function OptAtePairing(Pt:FpPoint;Qt:Fp3Point):Fp18Int;  // Compute Optimal Ate Pairings
              procedure SetLoopMode(Value:LoopPoweringMode);
              procedure SetPArams(const Value: StandardKSS18Curves);
            public
              SecurityLevel:String;    /// A symbolic identifier of the curve
              PerferedPoweringMode:LoopPoweringMode;
              property CoordinatesSystem:CoordinatesSystem read getCoord write SetCoord;
              property Beta:integer read getBeta write SetBeta;
              property Sigma:string read Getsigma write SetSigma;
              property Gamma:string read GetGamma Write SetGamma;
              property u :string read Getu write Setu;
              property A:string read getA write setA;
              property B:string read getB write setB;
              property Atw :string read GetAtw;
              property Btw :string read GetBtw;
              property P:string read GetP;
              property H:string read GetH;
              property N:string read getN;
              property R:String read getR;
              property Rtw:string read getRtw;
              property Htw:string read getHtw;
              property Tr:string read getTr;
              property Lp:string read getLp;
              property TwistMode:TTwistModel read GetTwm;
              Constructor Create(AOwner : TComponent); override;
              destructor Destroy;
              procedure SetStandardCurveParametres(inParametres:StandardKSS18Curves);
              procedure SetCustomCurveParametres(inParametres:KSS18CurvesParamsDefinition);
              procedure GenerateParamsTreeView(Tree:TTreeView);Override;
              function Paire(P:FpPoint;Q:Fp3Point):Fp18Int;
              function GetRandomG1Point:FpPoint;
              function GetRandomG2Point:Fp3Point;
              function GetDefautG1Generator:FpPoint;
              function GetDefautG2Generator:Fp3Point;
              function HashToG1Point(id:string):FpPoint;
              function HashToG2Point(id:String):Fp3Point;
            published
              property Fp18PoweringMode:FpPoweringMode read _FP18PoweringMode write _FP18PoweringMode;
              property PairingAlgorithm:KSS18CurvesPairingAlgos read _PairingAlgorithm write _PairingAlgorithm ;
              property LoopMode:LoopPoweringMode read _LoopMode write SetLoopMode;
              property Parametres:StandardKSS18Curves read params write Setparams;
            end;

const
  KSS18_192_1_Params:KSS18CurvesParamsDefinition=(SecurityLevel:'192bit';u:'0x15000000007004210';Beta:-2;sigma:'0+u*1+u^2*0';A:'0';B:'2';TwistMode:twDType;GenratorX:'-1';TwistGeneratorSeed:55);
  KSS18_Razvan_256_Params:KSS18CurvesParamsDefinition=(SecurityLevel:'256bit';u:'0x3FFFFFFFFFFFFFFFFFFFFFFFFFFF7FFFFFFFFFFFFC00010';Beta:-2;sigma:'0+u*1+u^2*0';A:'0';B:'2';TwistMode:twDType;GenratorX:'-1';TwistGeneratorSeed:55);
  KSS18_Razvan_192_Params:KSS18CurvesParamsDefinition=(SecurityLevel:'192bit';u:'-0x2000000000000083FFFFC0';Beta:-2;sigma:'0+u*1+u^2*0';A:'0';B:'2';TwistMode:twDType;GenratorX:'-1';TwistGeneratorSeed:55);
  KSS18_Razvan_128_Params:KSS18CurvesParamsDefinition=(SecurityLevel:'128bit';u:'0x1000003FFE02';Beta:-3;sigma:'0+u*1+u^2*1';A:'0';B:'3';TwistMode:twDType;GenratorX:'1';TwistGeneratorSeed:55);
  KSS18_192_3_Params:KSS18CurvesParamsDefinition=(SecurityLevel:'192bit';u:'0xFFF8800010000000';Beta:-2;sigma:'0+u*1+u^2*1';A:'0';B:'2';TwistMode:twDType;GenratorX:'1';TwistGeneratorSeed:55);

  KSS18ParamsList:array[0..4] of string=('KSS18_Razvan_128_Params','KSS18_192_1_Params','KSS18_Razvan_192_Params','KSS18_192_3_Params','KSS18_Razvan_256_Params');
  KSS18ImplementedPairingAlgos:array[0..0] of string=('Optimal Ate Pairing');

  procedure ComputeKSS18Parametres(Params: KSS18CurvesParamsDefinition; var Result:PtrCurveParams);
  procedure Register;

implementation



Procedure FastLoadParams(param:StandardKSS18Curves;var Result:PtrCurveParams);
var rho,tmp,e:Lint;
    i:integer;
begin
if Result=nil then begin
                   new(Result);
                   new(Result.LoopBin);
                   new(Result.LoopNaf);
                   new(Result.LoopTateNaf);
                   new(Result.LoopTateBin);
                   end;
with Result^ do begin
case param of
scKSS18at256:begin
               SecurityLevel:='256bit';
               Family:=cfKSS18;
               u:='0x3FFFFFFFFFFFFFFFFFFFFFFFFFFF7FFFFFFFFFFFFC00010';
               P:='0xC30C30C30C30C30C30C30C30C2FFFFFFFFFFFFFF9E7A00F3CF3CF424924924924929E79D173CF3CF50F3C4279FEA3CD3CF4543CF3CF2D292CE24810522524B1CDC69AF4DE192BC55ABA96F2A961'
               +'A1ABA9505DB67E44571D2B5CB0AC912097A2EEEBD216D12F22226F8D69BB39D454EF3CA7157A6E0D4BACF932282659CC90D4BCB0A2DFD43737BD70772E8F6D633412DE770BDFDCB6F0BBF257639D2B8'
               +'453EF0E51477A7E683E1D06BB43A93C5F16EAE01A282A2226FBDFE78D';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xBF112A8AD278E8DCEBD930835BC44AA2B49E3A36F34FEE166402FC717231426CF7CD0FABB55D4B61D08356343F64C2023D364C20D6F2F6FF6DB6D5D4B543EBC7A6432C0DBCB299BFE106BD1D17'
               +'7054A432B3F499B3E345687E16A761B7BE521A18C1B3CE693F011210B620BF1C6CC00256B69099006BB2761663ED54BDC4187148894D7A6C2BA820E3001';
               Tr:='0x24924924924924924924924924912492492492491B6DB92492492495B6DB6DB6DB6DEDB6CDB6DB6DB7B248B6DB7B6DB66DB6F6DB6DB6D8024ADB6D7FFFF7249912476DB70376D9B6DBA4925B6'
               +'DA92495B6DB4912249256DB6B6D24B701';
               Lp:= '0x9249249249249249249249249247FFFFFFFFFFFFF6DB70';
               n:='0xC30C30C30C30C30C30C30C30C2FFFFFFFFFFFFFF9E7A00F3CF3CF424924924924929E79D173CF3CF50F3C4279FEA3CD3CF4543CF3CF2D292CE24810522524B1CDC69AF4DE192BC55ABA96F2A96'
               +'1A1ABA9505DB67E44571D2B5CB0AC90FC0559CA5988F23EE5FD902668D896A78B305CF38BA7C1497B0288637B4CB8A2F122E705E2EC04667F8574B99BB321B6F57D1C03095507DA6C154E72576C789'
               +'2720C81574DD0A0C78C998AAB4D9A84A6A83BA1CDF595D346BB8EBB308D';
               htw:='0x97B425ED097B425ED097B425ECF425ED097B425E25ED365ED097B59097B425ED0991ED0385097B42AA5EA0E5F34ECF2D09DAA5ED09705E64C6CFCF096B343A0C0335EE3DAEB26CC963328FC9'
               +'28BE40DD91A05B61D9077F98650AAB23D21A5B94A16EBE1C8CE82E2EE201742380600ACF64C52EBDB9AE74EA76E3F6A4A2A9A702ABF94371F2191DF733463A5B3D4C2307520EC105BFC5FC4FAF1DEF'
               +'E3BBA461D3713EB9041544BE9AD0AF9C9A8AE02FB7D4ECA78A77D10B6F98CED3DCD721DA3A6EC30770C2A1597C479096A7DEC799F17C1A156665DD48525E3673416F4C55DC3740829EC95CD6934614'
               +'5ECE7D80B9457C7003E79DC8EA9F793781282A534E0A1454536804FF8E29E1E045E51D046F9A383CD6A7305FF2E2868D8000DB00EF2AAD098F9B24B806681A40BF26CA34F31AE22416ACDEF9D0CEA4'
               +'2A8DD537A9468D97971D8C1A7401A872CD84FD149EB81069E4515F0B3DE6A8FC52B0AB3724E3EA6341161B38DCF3B57FF4969F332ADCA7513E891D2CD7F968D6594CD669717304E27F6CB4EDB38317'
               +'D16BA616E87BFC70DB3DDECBEEB2735D1EF5B906DA93DDB5D90B';
               rtw:=r;
               A:=0;
               B:=2;
               H:='0x105555555555555555555555555513FFFFFFFFFFFDF555EC6AAAAAAAEC0000000000041554272AAAAABAFFF68EAC08D';
               FieldParam:=InitFieldParamsFp3(-2,p);
               New(TowerParam3);
               TowerParam3^.FieldParam:=FieldParam;
               TowerParam3^.Sigma.SetFieldParams(TowerParam3^.FieldParam);
               TowerParam3.Sigma.SetFormString('0+u*1+u^2*0');
               TowerParam3.pmod8:=p mod 8;
               BtwFp3.SetFieldParams(FieldParam);
               AtwFp3.SetFieldParams(FieldParam);
               TwistMode:=twDType;
               BtwFp3.c:='0xC30C30C30C30C30C30C30C30C2FFFFFFFFFFFFFF9E7A00F3CF3CF424924924924929E79D173CF3CF50F3C4279FEA3CD3CF4543CF3CF2D292CE24810522524B1CDC69AF4DE192BC55ABA96'+'F2A961A1ABA9505DB67E44571D2B5CB0AC912097A2EEEBD216D12F22226F8D69BB39D454EF3CA7157A6E0D4BACF932282659CC90D4BCB0A2DFD43737BD70772E8F6D633412DE770BDFDCB6F0BBF257639D2B8453EF0E51477A7E683E1D06BB43A93C5F16EAE01A282A2226FBDFE78C';
               BtwFp3.a:=0;
               BtwFp3.b:=0;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointY:='1';
               BasePointX:='-1';
               TwistBasePointX_Fp3.SetFieldParams(TowerParam3^.FieldParam);
               TwistBasePointy_Fp3.SetFieldParams(TowerParam3^.FieldParam);
               TwistBasePointX_Fp3.a:='0xA15469569B135600BD707B422FD1303C408A38104DD373ACDCB847C1971A67486488897785342D2231B922FA066198F4D3F3A2B7978BFD11683F94F2BA6CEC7C81B455E4D9C'+'C9E335DB347566FB2F25CD16BCEBD95927F602EC22EFD5F052CB3F683755444D2E9332966EA6C36D923656C1A7F62B468292A0E351DCB4A0C4C9E4ECFEF253D379DDCCB70207B1C1539CB594628643AA1E26BCB3F4892A3DD686493332CAD13E524D944D9278F133D29CA1D00B06D1001E9B3C91';
               TwistBasePointX_Fp3.b:='0xA15469569B135600BD707B422FD1303C408A38104DD373ACDCB847C1971A67486488897785342D2231B922FA066198F4D3F3A2B7978BFD11683F94F2BA6CEC7C81B455E4D9C'+'C9E335DB347566FB2F25CD16BCEBD95927F602EC22EFD5F052CB3F683755444D2E9332966EA6C36D923656C1A7F62B468292A0E351DCB4A0C4C9E4ECFEF253D379DDCCB70207B1C1539CB594628643AA1E26BCB3F4892A3DD686493332CAD13E524D944D9278F133D29CA1D00B06D1001E9B3C8C';
               TwistBasePointX_Fp3.c:='0xA15469569B135600BD707B422FD1303C408A38104DD373ACDCB847C1971A67486488897785342D2231B922FA066198F4D3F3A2B7978BFD11683F94F2BA6CEC7C81B455E4D9C'+'C9E335DB347566FB2F25CD16BCEBD95927F602EC22EFD5F052CB3F683755444D2E9332966EA6C36D923656C1A7F62B468292A0E351DCB4A0C4C9E4ECFEF253D379DDCCB70207B1C1539CB594628643AA1E26BCB3F4892A3DD686493332CAD13E524D944D9278F133D29CA1D00B06D1001E9B3C8C';
               TwistBasePointY_Fp3.a:='0x5B47FB5BA8809AF108871043CCFA298464C92053ACC3861F06BAAA55F6AF357A1F98BAC7E6D9774451095C940E7F06B9B16FA805A8108931466BED8B27EF70AF1F680B7E3F9'+'2BCD1BB4FB1E46F1936149547F3292E28B848C1C59A3C37BDB4F361C7B55FD0A4654EA3F05EA9040F6A576F4A7DBAE8B88B9C1A1FBB14A459153FFFF6C8327AC362D8250C87E161D269111EF84310174FE3DF28F4D003035DF1A8F3553422B59136B17EAD885140D657046EDEF407B004FD623B9';
               TwistBasePointY_Fp3.b:='0x8B31576B9E1845391613548D9FA5D938EE9CFC126CA3A6E44E33F12CBAD560B8441625B81365B3F0F5F0BD5DAF1C66EEFB0582A0B9D82DDD08F04675CCE1403FDDC53C5FE08'+'4B76B09F161850FAA834A36956FC2BD217603949E17DA5810C35159F99EB640FFD686BD4EBB2C7B97748DE3BB551DE2303DC34C5B1047324C2997540C8AEFC16360F93CA9F7EFB9664612FC52BBBDD8C5404EAEE8DE8CD4F7A675E30930B67EE7F7EABC9EAF05B9DDA834321F54A5AF14AF8AD64';
               TwistBasePointY_Fp3.c:='0x8A7F639E63E602A8C07BB88EA7E2B4A4EE72C9B755EF684D1D515FC566C517E6D2339455BCEB422D1F6D51D4FD2FA6CABBB5224A4B8559DC6ECE2639CAEDF094EFD84271E72'+'018ADFC0AC2FAF0EEBB3EF45050994D7914138D33E2509880E362549700FAB1AC220A5093282ACF2E1DA65DE70A37849B480E418DC086C21D23A8BBA2175C4D6C44A9528D6766A2F8BF57005184C62147D084174AA41292BB95FB24CF192AECA66294CBD3D00E29C2445D0AE202B827214DCB98D';

                TowerParam3^.FrobeniusP_Const[0].SetFieldParams(TowerParam3^.FieldParam);
                TowerParam3^.FrobeniusP_Const[0]:=TowerParam3^.Sigma.Pow((P-1)/6);
                for i:=1 to 4 do begin
                                 TowerParam3^.FrobeniusP_Const[i].SetFieldParams(TowerParam3^.FieldParam);
                                 TowerParam3^.FrobeniusP_Const[i]:=TowerParam3^.FrobeniusP_Const[i-1]*TowerParam3^.FrobeniusP_Const[0];
                                 end;
                FrobeniusMapConstX_Fp3.SetFieldParams (TowerParam3.FieldParam);
                case TwistMode of
                    twDType: begin
                             FrobeniusMapConstX_Fp3:=TowerParam3^.Sigma.Pow((p-1)/3);
                             FrobeniusMapConstY_Fp3:=TowerParam3^.Sigma.Pow((p-1)/2);
                             end;
                twMType:begin
                        FrobeniusMapConstX_Fp3:=(TowerParam3^.Sigma.Pow((p-1)/3)).Inverse;///
                        FrobeniusMapConstY_Fp3:=(TowerParam3^.Sigma.Pow((p-1)/2)).Inverse;   //////
                        end;
                end;
        end;
scKSS18at192_1:begin
               SecurityLevel:='192bit';
               Family:=cfKSS18;
               u:='0x15000000007004210';
               P:='0x6B5A6E1D11E5108CACCD56558145C710554EEFF9AF413DF697C385A6E6516BAE211B5FF0FAB2EBAA158830E019751E0E2C2DAC611D8E26EC036085796084958D';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0x3D0BF00007A1C6041585C6109046D58B977C7F66A8313FC6FF461832480226A0CF5B1714CCE143790ADF7EC23E03001';
               Tr:='0x6C870000090B955A2C485F55BBED1956747085A39DBB9817430F9CFEF1379701';
               Lp:= '0x3000000001000970';
               n:='0x6B5A6E1D11E5108CACCD56558145C710554EEFF9AF413DF697C385A6E6516BADB4945FF0F1A7564FE93FD18A5D8804B7B7BD26BD7FD28ED4C050E87A6F4CFE8D';
               htw:='0x4F2ABC4ED60F17998147E554920B6B9BF6BF26A0EDF5A5BBF8AEDFBF6D5DE46D2F4851C68629DF8D2F0C192CE4E2748545A90666DAD311FB007F3D2B4E1C96'
               +'553C16E72B2D8A05CA50C170FA8A8158ABD9E3F5A9256E5F12A8770D9323BA1CC6C91725CDCE0E804C715600B84F5EFDDC08E841463D45CC04876D9821BD94C75035'
               +'1C1C5A2DF701900BBCEE0FAE07BCD0B';
               rtw:=r;
               A:=0;
               B:=2;
               H:='0x1C230000012C2B10D93320905AF2328E8D';
               FieldParam:=InitFieldParamsFp3(-2,p);
               New(TowerParam3);
               TowerParam3^.FieldParam:=FieldParam;
               TowerParam3^.Sigma.SetFieldParams(TowerParam3^.FieldParam);
               TowerParam3.Sigma.SetFormString('0+u*1+u^2*0');
               TowerParam3.pmod8:=p mod 8;
               BtwFp3.SetFieldParams(FieldParam);
               AtwFp3.SetFieldParams(FieldParam);
               TwistMode:=twDType;
               BtwFp3.c:='0x6B5A6E1D11E5108CACCD56558145C710554EEFF9AF413DF697C385A6E6516BAE211B5FF0FAB2EBAA158830E019751E0E2C2DAC611D8E26EC036085796084958C';
               BtwFp3.a:=0;
               BtwFp3.b:=0;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointY:='1';
               BasePointX:='-1';
               TwistBasePointX_Fp3.SetFieldParams(TowerParam3^.FieldParam);
               TwistBasePointy_Fp3.SetFieldParams(TowerParam3^.FieldParam);
               TwistBasePointX_Fp3.a:='0x4E6B1F3273FE73D6707906C337D70AF31F0810FD57AEA298431BC81D27D4F0055E3A5066583EFA7A211C78181A03EB3709204771529A5D7C54E210B15CEA2E8D';
               TwistBasePointX_Fp3.b:='0x4E6B1F3273FE73D6707906C337D70AF31F0810FD57AEA298431BC81D27D4F0055E3A5066583EFA7A211C78181A03EB3709204771529A5D7C54E210B15CEA2E89';
               TwistBasePointX_Fp3.c:='0x4E6B1F3273FE73D6707906C337D70AF31F0810FD57AEA298431BC81D27D4F0055E3A5066583EFA7A211C78181A03EB3709204771529A5D7C54E210B15CEA2E89';
               TwistBasePointY_Fp3.a:='0x1768307F828312DD91CA9BD9FDA5609A4A429BC303312386FBBB049E895DCD5D558BD96F739D2FB0E4614B7C3E617DC8FEDE8C2DCAB75B0098C9993FD38141C4';
               TwistBasePointY_Fp3.b:='0x31AB40CA45B89F7E67C21F2F186D28D9BFC927EF0887503EE19C267E7118D167525341DE708B8DA756E4335D2093B88B03181DDA7618327B56F0BBF108DA01BA';
               TwistBasePointY_Fp3.c:='0x3FB78F8486B1BD20427F4E5AC98935D4FF1925F98FEB659DE17D7BB94B6DF4493422F5E0A5F600D869DFC084762C34CE6175A1AB4DF0CF405551F748DDA64C66';

                TowerParam3^.FrobeniusP_Const[0].SetFieldParams(TowerParam3^.FieldParam);
                TowerParam3^.FrobeniusP_Const[0]:=TowerParam3^.Sigma.Pow((P-1)/6);
                for i:=1 to 4 do begin
                                 TowerParam3^.FrobeniusP_Const[i].SetFieldParams(TowerParam3^.FieldParam);
                                 TowerParam3^.FrobeniusP_Const[i]:=TowerParam3^.FrobeniusP_Const[i-1]*TowerParam3^.FrobeniusP_Const[0];
                                 end;
                FrobeniusMapConstX_Fp3.SetFieldParams (TowerParam3.FieldParam);
                case TwistMode of
                    twDType: begin
                             FrobeniusMapConstX_Fp3:=TowerParam3^.Sigma.Pow((p-1)/3);
                             FrobeniusMapConstY_Fp3:=TowerParam3^.Sigma.Pow((p-1)/2);
                             end;
                twMType:begin
                        FrobeniusMapConstX_Fp3:=(TowerParam3^.Sigma.Pow((p-1)/3)).Inverse;///
                        FrobeniusMapConstY_Fp3:=(TowerParam3^.Sigma.Pow((p-1)/2)).Inverse;   //////
                        end;
                end;
        end;
scKSS18at192_2:begin
               SecurityLevel:='192bit';
               Family:=cfKSS18;
               u:='-0x2000000000000083FFFFC0';
               P:='0xC30C30C30C30DC30C2FFE18619CC8616FE1618C43A2DA5ACE3BF81206305E834BC180EE15B76CD748BD5DCFBF65CCC35F527C35E771FBD43616D4719941AAB02643C88002A0A2030CF17D0E587723453FF33EA0DD';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0x2FC44AA2B49E3ED5752B49E3A3A34F1CED6F3509FABC72F1782A87A4C0981692AD2959E9865E24242D3CE8B045AF49DF91E16BDFDD482F4CD86EC6C38C0001';
               Tr:='0x24924924924926EDB6DA49249257BA49166DB6DF03896D7C84926597136B1D7ADBA80DB6C892364B6DC01';
               Lp:= '0x4924924924924A5249240';
               n:='0xC30C30C30C30DC30C2FFE18619CC8616FE1618C43A2DA5ACE3BF81206305E834BC180EE15B76CD748BD5B869AD3839ECCE3A0C842DFB2AEBA72430ABDD3BA778F6C0036DC4730CC5B19CF53D79BB6BC1C8E87C4DD';
               htw:='0x25ED097B425EDB97B42085ED0AFCF25D55B909F972294A6BAAEA952A6510676F3ADE7523BFE03BDB53710A6A8FFE29C43D6D5128351F58F04A9A61218BB347F0E31CFED057D3BEBD0CAFC44911E'
               +'4FE6061C1BCB0DFE1CABD6C1018C4723C40F7515047AC02CCDFA4CA67012C0C71308CCC5D9EA0D1522928DD47C1360CFBB2F9EB7A20C9B02FA315F92EDAE32FD21BF0767C6781A72F9FCA09390A5A5A62'
               +'FAD603D3C5A297B789B509ECB30E4F53D5BE6EB28D866CBCA7497A05C3783969EB';
               rtw:=r;
               A:=0;
               B:=2;
               H:='0x4155555555555770555445CAAAAF025AA64A8EABC4DD';
               FieldParam:=InitFieldParamsFp3(-2,p);
               New(TowerParam3);
               TowerParam3^.FieldParam:=FieldParam;
               TowerParam3^.Sigma.SetFieldParams(TowerParam3^.FieldParam);
               TowerParam3.Sigma.SetFormString('0+u*1+u^2*0');
               TowerParam3.pmod8:=p mod 8;
               BtwFp3.SetFieldParams(FieldParam);
               AtwFp3.SetFieldParams(FieldParam);
               TwistMode:=twDType;
               BtwFp3.c:='0xC30C30C30C30DC30C2FFE18619CC8616FE1618C43A2DA5ACE3BF81206305E834BC180EE15B76CD748BD5DCFBF65CCC35F527C35E771FBD43616D4719941AAB02643C88002A0A2030CF17D0E587723453FF33EA0DC';
               BtwFp3.a:=0;
               BtwFp3.b:=0;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointY:='1';
               BasePointX:='-1';
               TwistBasePointX_Fp3.SetFieldParams(TowerParam3^.FieldParam);
               TwistBasePointy_Fp3.SetFieldParams(TowerParam3^.FieldParam);
               TwistBasePointX_Fp3.a:='0x7C22864A7199D843C8E818E1E4F8F2794292AA6E8A1918D0257C9DCEA13911F388841E0E1EC436618FDE7EE6C056656C9C67A3979EEAC35AD06FD6AA40412CA721BD1866E67CD4BED46C35174A3991617D96B1EBD';
               TwistBasePointX_Fp3.b:='0x7C22864A7199D843C8E818E1E4F8F2794292AA6E8A1918D0257C9DCEA13911F388841E0E1EC436618FDE7EE6C056656C9C67A3979EEAC35AD06FD6AA40412CA721BD1866E67CD4BED46C35174A3991617D96B1EBB';
               TwistBasePointX_Fp3.c:='0x7C22864A7199D843C8E818E1E4F8F2794292AA6E8A1918D0257C9DCEA13911F388841E0E1EC436618FDE7EE6C056656C9C67A3979EEAC35AD06FD6AA40412CA721BD1866E67CD4BED46C35174A3991617D96B1EBB';
               TwistBasePointY_Fp3.a:='0x41B1056995942D500A06F934B7273CCCEDBA9D319C9DA90BB5715485A338CF6E5659FD2D4E29E7335E2F9A2F85D67DDE8C058A95FF6AC9A4C47E0E47593B2C41FB97941CB620962070B9FAE4630D360C29B092ACF';
               TwistBasePointY_Fp3.b:='0x8D7846B69581222A6018DEA03C4D816F4EC618E5235F400DE03CDF521068AD9F002F41A6293B9F6D004425AF5596C1DA2B30EDE70E88F024CE07DC1E493B34B21F9647DAC733C8C6065072B3B4803098F6A401766';
               TwistBasePointY_Fp3.c:='0xAF760A7E1AF736CD188B190DB8A822C29B0A89C94B87EF77F95C1E37AD1273149115C3D496A40DB7463EC1A861792602D86DEE99FD1D43AA9828522C7A6FC648DC0250B795BB94030B4B640FF22F920F0AB6A5B7E';

                TowerParam3^.FrobeniusP_Const[0].SetFieldParams(TowerParam3^.FieldParam);
                TowerParam3^.FrobeniusP_Const[0]:=TowerParam3^.Sigma.Pow((P-1)/6);
                for i:=1 to 4 do begin
                                 TowerParam3^.FrobeniusP_Const[i].SetFieldParams(TowerParam3^.FieldParam);
                                 TowerParam3^.FrobeniusP_Const[i]:=TowerParam3^.FrobeniusP_Const[i-1]*TowerParam3^.FrobeniusP_Const[0];
                                 end;
                FrobeniusMapConstX_Fp3.SetFieldParams (TowerParam3.FieldParam);
                case TwistMode of
                    twDType: begin
                             FrobeniusMapConstX_Fp3:=TowerParam3^.Sigma.Pow((p-1)/3);
                             FrobeniusMapConstY_Fp3:=TowerParam3^.Sigma.Pow((p-1)/2);
                             end;
                twMType:begin
                        FrobeniusMapConstX_Fp3:=(TowerParam3^.Sigma.Pow((p-1)/3)).Inverse;///
                        FrobeniusMapConstY_Fp3:=(TowerParam3^.Sigma.Pow((p-1)/2)).Inverse;   //////
                        end;
                end;
        end;
scKSS18at192_3:begin
               SecurityLevel:='256bit';
               Family:=cfKSS18;
               u:='0xFFF8800010000000';
               P:='0xC2DE7E97B99E94E45AF42ED725E7A54AF17262FCF4EABB1E5CE0BFDC37FD978B1187740D742A4585A55B65AEFC72FC968D5EE18CB60B05401AEEDBC1AAAAB1D';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xBEEF96FC43E8FC6CC5844209890846FAF0455D20D96954701574E0F6A1ACD51A9EF6DB0E3000000000000000000001';
               Tr:='0x248E00303F0E21CA09F0E200303FFFBB6DB700000000000249136DB700000001';
               Lp:= '0x249136DB70000000';
               n:='0xC2DE7E97B99E94E45AF42ED725E7A54AF17262FCF4EABB1E5CE0BFDC37FD9788C8A77109834828E5064D45ABF87300DFB1EEE18CB60B051B89B80051AAAAB1D';
               htw:='0x976439D6974885E9C080BCC2317DB95C73D67292F410CE2AF0872800262ED4D0E336346A50D2B8F9BD2CF9F7DA0133832492C76F59CEF5FAD24FA26429AAD9'
               +'188B1F1281C43B6AF2D9EE8563459009556B2EC8F6E4E0956E75593556FD745199917E55B8003496115B779125AA237A9E2DB6DBD67AF4518972F2F1846CD296E8E3'
               +'20BB7CA9A00ADC5DF8BA9B30D4B6B';
               rtw:=r;
               A:=0;
               B:=2;
               H:='0x10546058EE1FF0B051B89B80051AAAAB1D';
               FieldParam:=InitFieldParamsFp3(-2,p);
               New(TowerParam3);
               TowerParam3^.FieldParam:=FieldParam;
               TowerParam3^.Sigma.SetFieldParams(TowerParam3^.FieldParam);
               TowerParam3.Sigma.SetFormString('0+u*1+u^2*1');
               TowerParam3.pmod8:=p mod 8;
               BtwFp3.SetFieldParams(FieldParam);
               AtwFp3.SetFieldParams(FieldParam);
               TwistMode:=twDType;
               BtwFp3.c:='1';
               BtwFp3.a:='2';
               BtwFp3.b:='0xC2DE7E97B99E94E45AF42ED725E7A54AF17262FCF4EABB1E5CE0BFDC37FD978B1187740D742A4585A55B65AEFC72FC968D5EE18CB60B05401AEEDBC1AAAAB1B';
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointY:='1';
               BasePointX:='-1';
               TwistBasePointX_Fp3.SetFieldParams(TowerParam3^.FieldParam);
               TwistBasePointy_Fp3.SetFieldParams(TowerParam3^.FieldParam);
               TwistBasePointX_Fp3.a:='0x557AFB98E62FC00CE5D753289A02CF7047D2BDE7BD69C6CE047802A92D5D73137A784E14CAF2068031A31F67B58CC7E941CB2DC8E563B844AC87E48BCEA2BDF';
               TwistBasePointX_Fp3.b:='0x557AFB98E62FC00CE5D753289A02CF7047D2BDE7BD69C6CE047802A92D5D73137A784E14CAF2068031A31F67B58CC7E941CB2DC8E563B844AC87E48BCEA2BDB';
               TwistBasePointX_Fp3.c:='0x557AFB98E62FC00CE5D753289A02CF7047D2BDE7BD69C6CE047802A92D5D73137A784E14CAF2068031A31F67B58CC7E941CB2DC8E563B844AC87E48BCEA2BDB';
               TwistBasePointY_Fp3.a:='0xBA448A0ECC04FE203C9546255F064FA6CEA30D9BCA59960D1A77D4A070E23E1BED296D59661EA64B66AD945B07B8322ED768ABD0BDFF24DCC32B753CE869EC1';
               TwistBasePointY_Fp3.b:='0xAB1244AAA09651767EA1BC1EFDFBA56E81C8CEE7761ADE344EAD803F7D9EC0F748DFF19FD4DBE017EEBD23BD89E0B7E9B4061409514AA619AC06A4E9977FBB';
               TwistBasePointY_Fp3.c:='0x3F1D72D15F321542047732D019516A53314CA87DBF623729F6E6764905A6D370474B871FCE976E452B765FA384BC97D5E7B7F6B82EDD1D0EE6983B54E651DDB';

                TowerParam3^.FrobeniusP_Const[0].SetFieldParams(TowerParam3^.FieldParam);
                TowerParam3^.FrobeniusP_Const[0]:=TowerParam3^.Sigma.Pow((P-1)/6);
                for i:=1 to 4 do begin
                                 TowerParam3^.FrobeniusP_Const[i].SetFieldParams(TowerParam3^.FieldParam);
                                 TowerParam3^.FrobeniusP_Const[i]:=TowerParam3^.FrobeniusP_Const[i-1]*TowerParam3^.FrobeniusP_Const[0];
                                 end;
                FrobeniusMapConstX_Fp3.SetFieldParams (TowerParam3.FieldParam);
                case TwistMode of
                    twDType: begin
                             FrobeniusMapConstX_Fp3:=TowerParam3^.Sigma.Pow((p-1)/3);
                             FrobeniusMapConstY_Fp3:=TowerParam3^.Sigma.Pow((p-1)/2);
                             end;
                twMType:begin
                        FrobeniusMapConstX_Fp3:=(TowerParam3^.Sigma.Pow((p-1)/3)).Inverse;///
                        FrobeniusMapConstY_Fp3:=(TowerParam3^.Sigma.Pow((p-1)/2)).Inverse;   //////
                        end;
                end;
        end;
scKSS18at128:begin
               SecurityLevel:='256bit';
               Family:=cfKSS18;
               u:='0x1000003FFE02';
               P:='0xC30C4923D19233F63B15EC76A0B1E3811876D305D6C7F094EBD7A3B4C236E522899C0A653BA5D7DE195C70F';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xBF113C73E06B6E8506050EB690BC396C6AC0BDBD55D5EB0CE15B42022EA458D9';
               Tr:='0x24924B6DA4B2484A27F5BDF69FD7EC6D3A47F89FF251';
               Lp:='0x249249B6D6E';
               n:='0xC30C4923D19233F63B15EC76A0B1E3811876D305D6C5A77034FD58903D9465C6AA320CE674D2335E8F5D4BF';
               htw:='0x97B45096665E76E7BA3CACEE6B0B13DEDD0E65B2C1C0847EADC0958A7A14C1A9DE3CEACBB8ADC41ACACE991'
               +'22E595FBED8495F73066DECCBA858181C216CE6A8D97DA25355C8773E8FFD205B45B657AE94BD52F365DF6CE90BE0E6B93C315207B6BFF';
               rtw:=r;
               A:=0;
               B:=3;
               H:='0x105555D7FBF4DFEFCFAADB57';
               FieldParam:=InitFieldParamsFp3(-3,p);
               New(TowerParam3);
               TowerParam3^.FieldParam:=FieldParam;
               TowerParam3^.Sigma.SetFieldParams(TowerParam3^.FieldParam);
               TowerParam3.Sigma.SetFormString('0+u*1+u^2*1');
               TowerParam3.pmod8:=p mod 8;
               BtwFp3.SetFieldParams(FieldParam);
               AtwFp3.SetFieldParams(FieldParam);
               TwistMode:=twDType;
               BtwFp3.c:='0x61862491E8C919FB1D8AF63B5058F1C08C3B6982EB63F84A75EBD1DA611B729144CE05329DD2EBEF0CAE388';
               BtwFp3.a:='0x61862491E8C919FB1D8AF63B5058F1C08C3B6982EB63F84A75EBD1DA611B729144CE05329DD2EBEF0CAE389';
               BtwFp3.b:='0x61862491E8C919FB1D8AF63B5058F1C08C3B6982EB63F84A75EBD1DA611B729144CE05329DD2EBEF0CAE386';
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointY:='1';
               BasePointX:='-1';
               TwistBasePointX_Fp3.SetFieldParams(TowerParam3^.FieldParam);
               TwistBasePointy_Fp3.SetFieldParams(TowerParam3^.FieldParam);
               TwistBasePointX_Fp3.a:='0x2594298976BD1616DFB589177EC7CDDDD8AFE279324FED6F2D59397F4284332366618F00AC9824031D1BD21';
               TwistBasePointX_Fp3.b:='0x2594298976BD1616DFB589177EC7CDDDD8AFE279324FED6F2D59397F4284332366618F00AC9824031D1BD20';
               TwistBasePointX_Fp3.c:='0x2594298976BD1616DFB589177EC7CDDDD8AFE279324FED6F2D59397F4284332366618F00AC9824031D1BD20';
               TwistBasePointY_Fp3.a:='0x697A1579D5D8685BEE0CD7DE50EFFC3EB81C368435B53686CE7A95EAB3AE25FEE13B3309A64E12DDA37BAF1';
               TwistBasePointY_Fp3.b:='0x48370B772C623C2CD390CF4BA62F777E869EC540A24FA4F68628DE75CD560A310B3E27DB1FC257802906F3C';
               TwistBasePointY_Fp3.c:='0x42B7D8E8B7028BC60C70697DFF3A1D61CDE8285E2B890F4B8F7B425AFE641ECB8037F0B75CC99C8E82E5284';

                TowerParam3^.FrobeniusP_Const[0].SetFieldParams(TowerParam3^.FieldParam);
                TowerParam3^.FrobeniusP_Const[0]:=TowerParam3^.Sigma.Pow((P-1)/6);
                for i:=1 to 4 do begin
                                 TowerParam3^.FrobeniusP_Const[i].SetFieldParams(TowerParam3^.FieldParam);
                                 TowerParam3^.FrobeniusP_Const[i]:=TowerParam3^.FrobeniusP_Const[i-1]*TowerParam3^.FrobeniusP_Const[0];
                                 end;
                FrobeniusMapConstX_Fp3.SetFieldParams (TowerParam3.FieldParam);
                case TwistMode of
                    twDType: begin
                             FrobeniusMapConstX_Fp3:=TowerParam3^.Sigma.Pow((p-1)/3);
                             FrobeniusMapConstY_Fp3:=TowerParam3^.Sigma.Pow((p-1)/2);
                             end;
                twMType:begin
                        FrobeniusMapConstX_Fp3:=(TowerParam3^.Sigma.Pow((p-1)/3)).Inverse;///
                        FrobeniusMapConstY_Fp3:=(TowerParam3^.Sigma.Pow((p-1)/2)).Inverse;   //////
                        end;
                end;
        end;
end;
end;
end;

procedure ComputeKSS18Parametres(Params: KSS18CurvesParamsDefinition; var Result:PtrCurveParams);
var u4,u6,u3,tmp,tmp1:LInt;
    i:integer;
    m2:HLInt;
    e:Lint;
    twtmp,twtmp1:Fp3Int;
    tau:array[0..4] of HLint;
    Genrator,Gtmp:Fp3Point;
begin
if Result=nil then begin
                   new(Result);
                   new(Result.LoopBin);
                   new(Result.LoopNaf);
                   new(Result.LoopTateNaf);
                   new(Result.LoopTateBin);
                   end;
with Result^ do begin
                SecurityLevel:=Params.SecurityLevel;
                Family:=cfKSS18;
                if (Pos('0x',Params.u)<>0)or(Pos('$',Params.u)<>0) then u:=_hex_To_LInt(Params.u)
                else u:=_str_To_LInt(Params.u);
                u4:=(u.Sqr).Sqr;
                u3:=u.Sqr*u;
                u6:=u4*u.Sqr;
                Tr:=(u4 + 16*u + 7)/7;
                r:=(u6 + 37*u3 + 343)/343;
                h:=(49*u.Sqr+245*u+343)/3;
                p:=h*r+tr-1;
                n:=h*r;
                if (not IsLIntPrime(p))or(not IsLIntPrime(r)) then
                    raise Exception.Create('Paramétre t invalide pour la construction de la courbe KSS...')
                else begin
                     p.Limit:=2*p.Data.i16[-2];
                     new(p.InverseFordivision);
                     p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
                     p.id:='Y';
                     Lp:=(u/7).Absolute;
                     end;
                A:=_Str_To_LInt(Params.A);
                B:=_Str_To_LInt(Params.B);
                FieldParam:=InitFieldParamsFp3(Params.Beta,p);
                //computing number of points on E(FP3)with cofactor
                Tau[0]:=tr.Sqr*tr-3*p*tr;
                Tau[1]:=p.Sqr*p;
                _Sqrt_HLint(((4*Tau[1]-Tau[0].Sqr)/3),Tau[2]);

                m2:=Tau[1]+1-((Tau[0]+3*Tau[2])/2);
                _VHCopy_LInt(m2,tmp1);
                Htw:=tmp1/r;    // Twist Cofactor
                {Htw:=u.Pow(18)+15*u.Pow(17)+96*u.Pow(16)+409*u.Pow(15)+1791*u.Pow(14)+7929*u.Pow(13)+
                 27539*u.Pow(12)+ 81660*u.Pow(11) + 256908*u.Pow(10) + 757927*u.Pow(9) + 1803684*u.Pow(8)+
                 4055484*u.Pow(7) + 9658007*u.Pow(6) + 19465362*u.Pow(5) + 30860595*u.Pow(4)+ 50075833*u.Pow(3)
                 + 82554234*u.Sqr + 88845918*u + 40301641;  }
                Rtw:=r;         ///  The true order is m2=r*htw , so we can choose a subgourp of order R (the cofactor is htw) (https://eprint.iacr.org/2005/133.pdf)

                New(TowerParam3);
                TowerParam3^.FieldParam:=FieldParam;
                TowerParam3^.Sigma.SetFieldParams(TowerParam3^.FieldParam);
                TowerParam3.Sigma.SetFormString(Params.Sigma);
                     ///   should test if beta and sigma are valide parametres
                Twtmp.SetFieldParams(FieldParam);

                // In Order to make this test , you have to incrase integers precision, but this slow down computation, so we consider that paramters
                // are already cheked befor input
                {_Pow_FP3(TowerParam3.Sigma,(Result.P.Sqr*P-1)/2,Twtmp);
                if (Twtmp.a=1)and(Twtmp.b=0) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe KSS...Sigma is a square');
                _Pow_FP3(Result.TowerParam3.Sigma,(Result.P.Sqr*P-1)/3,Twtmp);
                if (Twtmp.a=1)and(Twtmp.b=0) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe KSS...Sigma is a Cube');}
                    ///
                TowerParam3.pmod8:=p mod 8;
                BtwFp3.SetFieldParams(FieldParam);
                AtwFp3.SetFieldParams(FieldParam);
                TwistMode:=Params.TwistMode;
                if TwistMode=twMType then _Mul_FP_FP3(B,TowerParam3.Sigma,BtwFp3) else
                _Mul_FP_FP3(B,TowerParam3.Sigma.Inverse,BtwFp3);

                LoopBin^:=Lp.ToIntArray;
                LoopNaf^:=Lp.ToNafArray;
                LoopTateBin^:=n.ToIntArray;
                LoopTateNaf^:=n.ToNafArray;
                if (Pos('0x',Params.GenratorX)<>0)or(Pos('$',Params.GenratorX)<>0) then BasePointX:=_hex_To_LInt(Params.GenratorX)
                else BasePointX:=_Str_To_LInt(Params.GenratorX);
                ///   Constructing the base point of the generator for G1
                _Sqr_LInt(BasePointX,tmp);
                _Add_LInt(tmp,A,tmp);
                _Mul_LInt(tmp,BasePointX,tmp1);
                _Add_LInt(tmp1,B,tmp);
                _Mod_LInt(tmp,P,tmp);
                if tmp.IsASqrMod(P,FieldParam.MontgomeryData)
                                              then begin
                                                   ModSquareLInt(tmp,P, BasePointY,FieldParam.MontgomeryData,false);
                                                   if P-BasePointY<BasePointY then BasePointY:=P-BasePointY;
                                                   end
                else raise Exception.Create('Invalide parameter for the generator of the curve');
                ///   Constructing the base point for G2 Generator
                TwistBasePointX_Fp3.SetFieldParams(TowerParam3^.FieldParam);
                TwistBasePointy_Fp3.SetFieldParams(TowerParam3^.FieldParam);

                TwistBasePointX_Fp3.a.SetToRandom(TowerParam3^.FieldParam.p,params.TwistGeneratorSeed);
                TwistBasePointX_Fp3.b.SetToRandom(TowerParam3^.FieldParam.p,params.TwistGeneratorSeed);
                TwistBasePointX_Fp3.c.SetToRandom(TowerParam3^.FieldParam.p,params.TwistGeneratorSeed);

                repeat
                _Inc_LInt(TwistBasePointX_Fp3.a,1);
                _Sqr_FP3(TwistBasePointX_Fp3,twtmp);
                _Add_FP3(twtmp,AtwFp3,twtmp);
                _Mul_FP3(twtmp,TwistBasePointX_Fp3,twtmp1);
                _Add_FP3(twtmp1,BtwFp3,twtmp);
                until twtmp.IsASquare;
                _Sqrt_FP3(twtmp,TwistBasePointy_Fp3);
                _Neg_FP3(TwistBasePointy_Fp3,twtmp1);
                if _Compare_FP3(TwistBasePointy_Fp3,twtmp1)=1 then TwistBasePointy_Fp3:=twtmp1;



                TowerParam3^.FrobeniusP_Const[0].SetFieldParams(TowerParam3^.FieldParam);
                TowerParam3^.FrobeniusP_Const[0]:=TowerParam3^.Sigma.Pow((P-1)/6);
                for i:=1 to 4 do begin
                                 TowerParam3^.FrobeniusP_Const[i].SetFieldParams(TowerParam3^.FieldParam);
                                 TowerParam3^.FrobeniusP_Const[i]:=TowerParam3^.FrobeniusP_Const[i-1]*TowerParam3^.FrobeniusP_Const[0];
                                 end;
                FrobeniusMapConstX_Fp3.SetFieldParams (TowerParam3.FieldParam);
                case TwistMode of
                    twDType: begin
                             FrobeniusMapConstX_Fp3:=TowerParam3^.Sigma.Pow((p-1)/3);
                             FrobeniusMapConstY_Fp3:=TowerParam3^.Sigma.Pow((p-1)/2);
                             end;
                twMType:begin
                        FrobeniusMapConstX_Fp3:=(TowerParam3^.Sigma.Pow((p-1)/3)).Inverse;///
                        FrobeniusMapConstY_Fp3:=(TowerParam3^.Sigma.Pow((p-1)/2)).Inverse;   //////
                        end;
                end;
                Genrator.SetCurveParams(Result);
                Gtmp.SetCurveParams(Result);
                Gtmp.X:=TwistBasePointX_Fp3;
                Gtmp.Y:=TwistBasePointy_Fp3;
                Gtmp.Z.a:='1';
                Gtmp.Z.b:='0';
                Gtmp.Z.c:='0';
                Gtmp.Infinity:=false;
                _Mul_Fp_Fp3Point(htw,Gtmp,Genrator);
                TwistBasePointX_Fp3:=Gtmp.X;
                TwistBasePointy_Fp3:=Gtmp.Y;
                end;
end;




{ TKSS18Curve }

procedure TKSS18CurvePairing.GenerateParamsTreeView(Tree:TTreeView);
var s:string;
    r:real;
    G:FpPoint;Gtw:Fp3Point;
    tmp:TTreeNode;
begin
Tree.Items.Clear;
tmp:=Tree.Items.Add(nil,'Curve');
Tree.Items.AddChild(tmp,'Family : KSS18');
if CurveParams.A<0 then  s:='y^2=x^3-'+CurveParams.B.Absolute.ToDecimalString
else s:='y^2=x^3+'+CurveParams.B.ToDecimalString;
Tree.Items.AddChild(tmp,'Equation : '+s);
tmp:=Tree.Items.Add(nil,'Initial Parameter x0');
Tree.Items.AddChild(tmp,'x0= '+curveparams.u.ToHexString );
Tree.Items.AddChild(tmp,'Hamming Weight (Bin) ='+IntToStr(HammingWeight(CurveParams.u,False)));
Tree.Items.AddChild(tmp,'Hamming Weight (Neg) ='+IntToStr(HammingWeight(CurveParams.u,True)));
tmp:=Tree.Items.Add(nil,'Field Fp');
Tree.Items.AddChild(tmp,'Prime P = '+CurveParams.P.ToHexString );
Tree.Items.AddChild(tmp,'Size of Fp : '+inttostr(CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Field E(Fp)');
Tree.Items.AddChild(tmp,'Order (#E(Fp)) N = '+CurveParams.N.ToHexString );
Tree.Items.AddChild(tmp,'Frobenius Trace Tr : '+CurveParams.Tr.ToHexString);
tmp:=Tree.Items.Add(nil,'Base Fields G1 : E(Fp)[R]');
Tree.Items.AddChild(tmp,'Order of G1 R = '+CurveParams.R.ToHexString);
Tree.Items.AddChild(tmp,'Size of G1 : '+inttostr(CurveParams.R.BitLength)+' bit');
Tree.Items.AddChild(tmp,'Cofactor of G1 H = '+CurveParams.H.ToHexString);
tmp:=Tree.Items.Add(nil,'Twist');
Tree.Items.AddChild(tmp,'Equation of the Twist : '+'y^2=x^3+'+Curveparams.BtwFp4.toHexString);
Tree.Items.AddChild(tmp,'Degree of the Twist : Sextic');
if Curveparams.TwistMode=twMType then S:='M-Type' else S:='D-Type';
Tree.Items.AddChild(tmp,'Type of the Twist : '+S);
tmp:=Tree.Items.Add(nil,'Field G2 (Twist) :E''(Fp3)[R]');
Tree.Items.AddChild(tmp,'Order of G2 R = '+CurveParams.Rtw.ToHexString);
Tree.Items.AddChild(tmp,'Size of G2 : '+inttostr(CurveParams.Rtw.BitLength)+' bit');
Tree.Items.AddChild(tmp,'Cofactor of G2 Htw = '+CurveParams.Htw.ToHexString);
tmp:=Tree.Items.Add(nil,'Targted Field (GT):');
Tree.Items.AddChild(tmp,'Extension Field : Fp18');
Tree.Items.AddChild(tmp,'Size of GT : '+inttostr(18*CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Security Level');
Tree.Items.AddChild(tmp,'Security in G1 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
Tree.Items.AddChild(tmp,'Security in G2 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
r:=1.92* exp((1/3)*ln(18*CurveParams.P.BitLength))*sqr(exp((1/3)*(ln(ln(18*CurveParams.P.BitLength)/ln(2)))));
Tree.Items.AddChild(tmp,'Security in GT (GNFS) : '+inttostr(Round(r))+' bit').ImageIndex:=2;
tmp:=Tree.Items.Add(nil,'Implemented Pairings');
Tree.Items.AddChild(tmp,'Otp-Ate');
tmp:=Tree.Items.Add(nil,'Tower Construction');
if Curveparams.FieldParam.Beta<0 then Tree.Items.AddChild(tmp,'Fp3<u>=ExstensionField<u,|u^3'+inttostr(Curveparams.FieldParam.Beta)+'>')
else Tree.Items.AddChild(tmp,'Fp3<u>=ExstensionField<u,|u^3+'+inttostr(Curveparams.FieldParam.Beta)+'>');
Tree.Items.AddChild(tmp,'Fp9<v>=ExstensionField<v,|v^3-('+Curveparams.TowerParam3.Sigma.toHexString+')>');
Tree.Items.AddChild(tmp,'Fp18<w>=ExstensionField<w,|w^3-v>');
tmp:=Tree.Items.Add(nil,' = '+floattostrf(CurveParams.P.BitLength/CurveParams.R.BitLength,ffGeneral, 2, 4 ));
tmp:=Tree.Items.Add(nil,'Default G1 Generator');
G.SetCurveParams (CurveParams);
G.SetToDefaultGenerator;
Gtw.SetCurveParams (CurveParams);
Gtw.SetToDefaultGenerator;
Tree.Items.AddChild(tmp,G.toHexString);
tmp:=Tree.Items.Add(nil,'Default G2 Generator');
Tree.Items.AddChild(tmp,Gtw.toHexString);
Tree.FullExpand;
Tree.Selected:=Tree.Items[0];
//Tree.SetFocus;
end;



{*******************************************************************************}
constructor TKSS18CurvePairing.Create;
begin
inherited create(Aowner);
CurveParams:=nil;
SetStandardCurveParametres(scKSS18at256);
CoordinatesSystem:=csProjective;
Fp18PoweringMode:=pmKarbina;
SetLoopMode(lpmAuto);
PairingAlgorithm:=KSS18pOptAte;
end;

{*******************************************************************************}
destructor TKSS18CurvePairing.Destroy;
begin
Dispose(CurveParams);
inherited;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getA: String;
begin
Result:=_A.ToHexString;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getAtw: String;
begin
Result:=CurveParams^.AtwFp3.toHexString;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getB: String;
begin
Result:=_B.ToHexString;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getBeta:integer;
begin
Result:=CurveParams^.FieldParam.Beta;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getBtw: String;
begin
Result:=CurveParams^.BtwFp3.toHexString;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getCoord: CoordinatesSystem;
begin
Result:=CurveParams^.CoordSys;
end;

function TKSS18CurvePairing.GetDefautG1Generator: FpPoint;
begin
Result.SetToDefaultGenerator;
end;

function TKSS18CurvePairing.GetDefautG2Generator: Fp3Point;
begin
Result.SetToDefaultGenerator;
end;

function TKSS18CurvePairing.getGamma: String;
begin
if  TwistMode=twDType then Result:=CurveParams^.TowerParam2.Gamma.toHexString
else Result:=CurveParams^.TowerParam2.Gamma.Inverse.toHexString
end;

{*******************************************************************************}
function TKSS18CurvePairing.getH: string;
begin
Result:=CurveParams^.H.toHexString;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getHtw: string;
begin
Result:=CurveParams^.Htw.toHexString;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getLp: string;
begin
Result:=CurveParams^.Lp.toHexString;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getN: string;
begin
Result:=CurveParams^.N.toHexString;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getP: string;
begin
Result:=CurveParams^.P.toHexString;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getR: string;
begin
Result:=CurveParams^.R.toHexString;
end;

function TKSS18CurvePairing.GetRandomG1Point: FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

function TKSS18CurvePairing.GetRandomG2Point: Fp3Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getRtw: String;
begin
Result:=CurveParams^.Rtw.toHexString;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getSigma: String;
begin
Result:=CurveParams^.TowerParam3.Sigma.toHexString
end;

{*******************************************************************************}
function TKSS18CurvePairing.getTr: string;
begin
Result:=CurveParams^.Tr.toHexString;
end;

{*******************************************************************************}
function TKSS18CurvePairing.GetTwm: TTwistModel;
begin
Result:=CurveParams.TwistMode;
end;

{*******************************************************************************}
function TKSS18CurvePairing.getu: String;
begin
Result:=_u.ToHexString;
end;

{*******************************************************************************}
function TKSS18CurvePairing.HashToG1Point(id: string): FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

{*******************************************************************************}
function TKSS18CurvePairing.HashToG2Point(id: String): Fp3Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

{*******************************************************************************}
procedure TKSS18CurvePairing.RecomputeParametres;
var tmp:KSS18CurvesParamsDefinition;
begin
tmp.u:=_u.toHexString;
tmp.Beta:=beta;
tmp.Sigma:=Sigma;
tmp.Gamma:=Gamma;
tmp.A:=_A.toHexString;
tmp.B:=_B.toHexString;
tmp.SecurityLevel:=SecurityLevel;
ComputeKSS18Parametres(tmp,CurveParams);
end;

{*******************************************************************************}
procedure TKSS18CurvePairing.SetStandardCurveParametres(inParametres:StandardKSS18Curves);
begin
case inParametres of
scKSS18at256: begin
              FastLoadParams(scKSS18at256,CurveParams);
              _A:=KSS18_Razvan_256_Params.A;
              _B:=_Str_To_LInt(KSS18_Razvan_256_Params.B);
              _u:=_Hex_To_LInt(KSS18_Razvan_256_Params.u);
              end;
scKSS18at192_1: begin
              FastLoadParams(scKSS18at192_1,CurveParams);
              _A:=KSS18_192_1_Params.A;
              _B:=_Str_To_LInt(KSS18_192_1_Params.B);
              _u:=_Hex_To_LInt(KSS18_192_1_Params.u);
              end;
scKSS18at192_2: begin
              FastLoadParams(scKSS18at192_2,CurveParams);
              _A:=KSS18_Razvan_192_Params.A;
              _B:=_Str_To_LInt(KSS18_Razvan_192_Params.B);
              _u:=_Hex_To_LInt(KSS18_Razvan_192_Params.u);
              end;
scKSS18at192_3: begin
              FastLoadParams(scKSS18at192_3,CurveParams);
              _A:=KSS18_192_3_Params.A;
              _B:=_Str_To_LInt(KSS18_192_3_Params.B);
              _u:=_Hex_To_LInt(KSS18_192_3_Params.u);
              end;
scKSS18at128:begin
              FastLoadParams(scKSS18at128,CurveParams);
              _A:=KSS18_Razvan_128_Params.A;
              _B:=_Str_To_LInt(KSS18_Razvan_128_Params.B);
              _u:=_Hex_To_LInt(KSS18_Razvan_128_Params.u);
              end;

end;
if HammingWeight(Lp,true)<HammingWeight(lp,false) then PerferedPoweringMode:=lpmNaf
else PerferedPoweringMode:=lpmBinary;
end;

{*******************************************************************************}
procedure TKSS18CurvePairing.SetA(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS18CurvePairing.SetB(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS18CurvePairing.SetBeta(Value: Integer);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS18CurvePairing.setCoord(value: CoordinatesSystem);
begin
CurveParams^.CoordSys:=value;
end;

{*******************************************************************************}
procedure TKSS18CurvePairing.SetCustomCurveParametres(inParametres: KSS18CurvesParamsDefinition);
begin
ComputeKSS18Parametres(inParametres,CurveParams);
end;

procedure TKSS18CurvePairing.SetGamma(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS18CurvePairing.SetSigma(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS18CurvePairing.Setu(Value: String);
begin
RecomputeParametres;
end;

function TKSS18CurvePairing.Paire(P: FpPoint; Q: Fp3Point): Fp18Int;
begin
if (P.Infinity)or(Q.Infinity) then raise Exception.Create('Can''t compute pairing for infinit points ....');
case PairingAlgorithm of
  KSS18pOptAte:Result:=OptAtePairing(P,Q);
end;
end;


{******************** Set the loop mode of the pairing :Binary / Negatve Representation*************}
procedure TKSS18CurvePairing.SetLoopMode(Value: LoopPoweringMode);
begin
_LoopMode:=Value;
if _LoopMode=lpmAuto then _LoopMode:=PerferedPoweringMode;
case _LoopMode of
lpmBinary:begin
          case PairingAlgorithm of
          KSS18pOptAte:Loop:=CurveParams.LoopBin;
          end;

          end;
lpmNaf:begin
       case PairingAlgorithm of
       KSS18pOptAte:Loop:=CurveParams.LoopNaf;
       end;
       end;
end;
end;

procedure TKSS18CurvePairing.SetPArams(const Value: StandardKSS18Curves);
var pr:TComponent;
    i:integer;
begin
SetStandardCurveParametres(Value);
pr:=Self.Owner;
for i:=0 to pr.ComponentCount-1 do begin
                                   if (pr.Components[i] is TTreeView1)and (TTreeView1(pr.Components[i]).Curve = self)
                                   then GenerateParamsTreeView(TTreeView(pr.Components[i]));
                                   end;
params:=value;
end;

{******************** Final Exponentiation Step :Compute f^((p^18-1)/r) ***********************}
Procedure TKSS18CurvePairing.FinalPowerOptimalAteKSS18(ff:Fp18Int;u:LInt;PoweringMode:FpPoweringMode;var Result:Fp18Int); // Compute f^((p^18-1)/r)
var xA,xB,t0,t1,t2,t3,t4,t5,t6,t7,tmp,f:Fp18Int;
    t:array[1..7]of Fp18Int;
begin
///***** Soft Exponentiation ***///
_Conjugate_FP18(ff,t[1]);
_Inv_FP18(ff,t[2]);
_Mul_FP18(t[1],t[2],f);
_Pow_FP18_P_i(f,3, t[3]);
_Mul_FP18(t[3],f,f);
///***** Hard Exponentiation ***///
_Pow_FP18(f,u,t[1],PoweringMode);
_Pow_FP18(t[1],u,t[2],PoweringMode);
_Pow_FP18(t[2],u,t[3],PoweringMode);
_Pow_FP18(t[3],u,t[4],PoweringMode);
_Pow_FP18(t[4],u,t[5],PoweringMode);
_Pow_FP18(t[5],u,t[6],PoweringMode);
_Pow_FP18(t[6],u,t[7],PoweringMode);
_Conjugate_FP18(t[1],tmp);
_Pow_FP18_P_i(tmp,2,xA);
_Conjugate_FP18(f,tmp);
_Pow_FP18_P_i(tmp,2,xB);
_Mul_FP18(xA,xB,t0);
_Conjugate_FP18(t[2],tmp);
_Pow_FP18_P_i(tmp,2,xB);
_Mul_FP18(t0,xB,t1);
_Mul_FP18(t0,t0,t0);
_Conjugate_FP18(f,tmp);
_Pow_FP18_P_i(tmp,2,xB);
_Mul_FP18(t0,xB,t0);
_Pow_FP18_P_i(t[1],1,xB);
_Mul_FP18(t0,xB,t0);
_Pow_FP18_P_i(t[2],5,xA);
_Pow_FP18_P_i(t[4],4,tmp);
_Mul_FP18(tmp,xA,xA);
_Conjugate_FP18(t[5],tmp);
_Pow_FP18_P_i(tmp,2,tmp);
_Mul_FP18(tmp,xA,xA);
_Mul_FP18(xA,xB,t5);
_Mul_FP18(t0,t0,t0);
_Mul_FP18(t0,t1,t3);
_Pow_FP18_P_i(t[1],5,xA);
_Conjugate_FP18(t[4],tmp);
_Pow_FP18_P_i(tmp,2,tmp);
_Mul_FP18(tmp,xA,xA);
_Pow_FP18_P_i(t[2],1,xB);
_Mul_FP18(xA,xB,t1);
xA:=xB;

_Sqr_FP18(xA,t0);
_Pow_FP18_P_i(t[2],4,xB);
_Mul_FP18(t0,xB,t0);
_Pow_FP18_P_i(t[1],4,xB);
_Mul_FP18(t3,xB,t2);
_Conjugate_FP18(t[1],tmp);
_Pow_FP18_P_i(tmp,2,xB);
_Mul_FP18(t3,xB,t4);
_Mul_FP18(t2,t2,t2);
_Conjugate_FP18(t[2],tmp);
_Pow_FP18_P_i(tmp,3,xB);
_Mul_FP18(t0,xB,t3);
_Conjugate_FP18(t[2],xB);
_Mul_FP18(t3,xB,t0);
_Mul_FP18(t3,t4,t4);
_Conjugate_FP18(t[3],tmp);
_Pow_FP18_P_i(tmp,3,xB);
_Mul_FP18(t0,xB,t0);
_Mul_FP18(t0,t2,t3);
_Pow_FP18_P_i(f,5,xB);
_Conjugate_FP18(t[3],tmp);
_Pow_FP18_P_i(tmp,2,tmp);
_Mul_FP18(tmp,xB,xB);
_Mul_FP18(t3,xB,t2);
_Mul_FP18(t3,t5,t3);
_Mul_FP18(t3,t2,t5);
_Conjugate_FP18(t[3],xB);
_Mul_FP18(t2,xB,t2);
_Conjugate_FP18(t[6],tmp);
_Pow_FP18_P_i(tmp,3,xA);
_Mul_FP18(xA,xB,t3);
_Mul_FP18(t2,t2,t2);
_Mul_FP18(t2,t4,t4);
_Pow_FP18_P_i(t[3],1,xB);
_Mul_FP18(t1,xB,t2);
xA:=xB;
_Conjugate_FP18(t[2],tmp);
_Pow_FP18_P_i(tmp,3,xB);
_Mul_FP18(xA,xB,t1);
_Mul_FP18(t2,t4,t6);
_Pow_FP18_P_i(T[4],1,xB);
_Mul_FP18(t2,xB,t4);
_Pow_FP18_P_i(t[3],4,xB);
_Mul_FP18(t6,xB,t2);
_Pow_FP18_P_i(t[5],4,xB);
_Conjugate_FP18(t[5],tmp);
_Pow_FP18_P_i(tmp,3,tmp);
_Mul_FP18(tmp,xB,xB);
_Mul_FP18(t6,xB,t7);
_Mul_FP18(t2,t4,t4);
_Pow_FP18_P_i(t[6],1,xB);
_Mul_FP18(t2,xB,t2);
_Mul_FP18(t4,t4,t4);
_Mul_FP18(t4,t5,t4);
_Conjugate_FP18(t[4],xA);
_Conjugate_FP18(t[4],tmp);
_Pow_FP18_P_i(tmp,3,xB);
_Mul_FP18(xA,xB,t5);
_Mul_FP18(t3,xB,t3);
_Pow_FP18_P_i(t[5],1,xA);
xB:=xA;
_Mul_FP18(xB,xA,t6);
_Mul_FP18(t6,t7,t7);
_Pow_FP18_P_i(f,3,xB);
_Mul_FP18(t5,xB,t6);
_Mul_FP18(t6,t4,t4);
_Conjugate_FP18(t[7],tmp);
_Pow_FP18_P_i(tmp,3,xB);
_Mul_FP18(t6,xB,t6);
_Mul_FP18(t4,t0,t0);
_Pow_FP18_P_i(t[6],4,xB);
_Mul_FP18(t4,xB,t4);
_Mul_FP18(t0,t0,t0);
_Conjugate_FP18(t[5],xB);
_Mul_FP18(t0,xB,t0);
_Mul_FP18(t7,t1,t1);
_Mul_FP18(t4,t7,t4);
_Mul_FP18(t1,t1,t1);
_Mul_FP18(t1,t2,t2);
_Mul_FP18(t0,t3,t1);
_Conjugate_FP18(t[3],tmp);
_Pow_FP18_P_i(tmp,3,xB);
_Mul_FP18(t1,xB,t0);
_Mul_FP18(t1,t6,t1);
_Mul_FP18(t0,t0,t0);
_Mul_FP18(t0,t5,t0);
_Conjugate_FP18(t[6],xB);
_Mul_FP18(t2,xB,t2);
_Mul_FP18(t2,t2,t2);
_Mul_FP18(t2,t4,t2);
_Mul_FP18(t0,t0,t0);
_Mul_FP18(t0,t3,t0);
_Mul_FP18(t2,t1,t1);
_Mul_FP18(t1,t0,t0);
_Mul_FP18(t1,xB,t1);
_Mul_FP18(t0,t0,t0);
_Mul_FP18(t0,t2,t0);
_Conjugate_FP18(t[7],tmp);
_Mul_FP18(tmp,f,xB);
_Mul_FP18(t0,xB,t0);
_Mul_FP18(t1,xB,t1);
_Mul_FP18(t0,t0,t0);
_Mul_FP18(t0,t1,Result);
end;

{******************** Compute the Optimal Ate Pairing ***********************}
Function TKSS18CurvePairing.OptAtePairing(Pt:FpPoint;Qt:Fp3Point):Fp18Int;
var T,Q1,Q2,mQt:Fp3Point;
    i:integer;
    f,f2,f3:Fp18Int;
begin
f.SetTowerParams(CurveParams.TowerParam3);
f.SetToOne;
T:=Qt;
_Neg_Fp3_Point(Qt,mQt);
T.SetPairingPointCoordinates(Pt.X,Pt.Y);
Q1.SetPairingPointCoordinates(Pt.X,Pt.Y);
Q2.SetPairingPointCoordinates(Pt.X,Pt.Y);
T.ComputeLigneValue:=true;
T.ComputeLigneAtFp6:=false;
for i:=Length(Loop^)-2 Downto 0 do begin
                        case CoordinatesSystem of
                          csAffine:_Double_Affine_Fp3_Point(T,T,true);
                          csJacobian:_Double_Jacobian_Fp3_Point(T,T,True);
                          csProjective:_Double_projective_Fp3_Point(T,T,True);
                        end;
                        _Sqr_FP18(f,f);
                        _Sparse_Mul_FP18(f,T.LineAtP,f, TwistMode);
                        if  Loop^[i]=1 then begin
                                         case CoordinatesSystem  of
                                         csAffine:_Add_Affine_Fp3_Point(Qt,T,T,True);
                                         csJacobian:_Add_Jacobian_Fp3_Point(Qt,T,T,True);
                                         csProjective:_Add_Projective_Fp3_Point(Qt,T,T,True);
                                         end;
                                         _Sparse_Mul_FP18(f,T.LineAtP,f,TwistMode);
                                         end
                        else if (Loop^[i]=-1) then begin
                                                case CoordinatesSystem of
                                                csAffine:_Add_Affine_Fp3_Point(mQt,T,T,True);
                                                csJacobian:_Add_Jacobian_Fp3_Point(mQt,T,T,True);
                                                csProjective:_Add_Projective_Fp3_Point(mQt,T,T,True);
                                                end;
                                                _Sparse_Mul_FP18(f,T.LineAtP,f,TwistMode);
                                                end;
                        end;

_Sqr_FP18(f,f2);
case CoordinatesSystem of
  csAffine:begin
           Q1:=T;
           _Double_Affine_Fp3_Point(T,T,true);
           end;
  csJacobian:begin
             _Jacobian_To_Affine_Fp3Point(T,Q1);
             _Double_Jacobian_Fp3_Point(T,T,True);
             end;
  csProjective:begin
               _Projective_To_Affine_Fp3Point(T,Q1);
               _Double_Projective_Fp3_Point(T,T,True);
               end;
end;
_Mul_FP18(f2,T.LineAtP,f2);
_Mul_FP18(f,f2,f3);
case CoordinatesSystem of
  csAffine:begin
           Q2:=T;
           _Add_Affine_Fp3_Point(Q1,T,T,True);
           end;
  csJacobian:begin
             _Jacobian_To_Affine_Fp3Point(T,Q2);
             _Add_Jacobian_Fp3_Point(Q1,T,T,True);
             end;
  csProjective:begin
               _Projective_To_Affine_Fp3Point(T,Q2);
               _Add_Projective_Fp3_Point(Q1,T,T,True);
               end;
end;
_Mul_FP18(f3,T.LineAtP,f3);
_Pow_FP18_P_i(f3,6,f3);
_Mul_FP18(f3,f2,f);
_Frobenius_Map_i(T,6,T);
case CoordinatesSystem of
  csAffine:_Add_Affine_Fp3_Point(Q2,T,T,True);
  csJacobian:_Add_Jacobian_Fp3_Point(Q2,T,T,True);
  csProjective:_Add_Projective_Fp3_Point(Q2,T,T,True);
end;
_Mul_FP18(f,T.LineAtP,f);
FinalPowerOptimalAteKSS18(f,CurveParams.u,Fp18PoweringMode,Result);
end;

procedure Register;
begin
  RegisterComponents('Pairings', [TKSS18CurvePairing]);
end;

end.
