unit FP9Arithmetic;

interface

uses System.SysUtils,VCL.dialogs, LargeIntegers, Fp3Arithmetic;


{*******************************************************************************************
   Développer par Faraoun Kamel Mohamed
   Université Djilali Liabes -Sidi Bel Abbes - Algérie
   kamel_mh@yahoo.fr

   Arithmetic computation over Fp9. Fp9 is the ninth extention of Fp (Tower extention of order 3 for FP3)
   with respect to the irriductible polynômial V^3-Sigma (V^3=Sigma). Elements are in the form a+b*V+c*V^2
   where a,b and c are from Fp3
********************************************************************************************}

type

  PtrTowerParams18=^TowerParams18;
  TowerParams18=record
              FieldParam:PtrFieldParams;    // Parmetres of the Base Field
              Sigma:Fp3Int;                 // Parametre of the Cubic Extention   (Fp->Fp3)
              pmod8:Lint;
              FrobeniusP_Const:array[0..4]of Fp3Int;  // Frobenius Constants for PreComputation of f^p, f^p^2, f^p^3 (KSS Curves)
              end;

  Fp9Int=record
  a,b,c:Fp3Int;
  Tower:PtrTowerParams18;
  public
    class operator Add(const Left, Right: Fp9Int): Fp9Int;
    class operator Subtract(const Left, Right: Fp9Int): Fp9Int;
    class operator Multiply(const Left, Right: Fp9Int): Fp9Int;
    class operator Multiply(const Left:Fp3Int; Right: Fp9Int): Fp9Int;
    class operator Negative(const Value: Fp9Int): Fp9Int;
    class operator Equal(const Left, Right: Fp9Int): Boolean;
    class operator NotEqual(const Left, Right: Fp9Int): Boolean;
    function ToHexString: string;
    function ToDecimalString: string;
    function Inverse:Fp9Int;
    function Pow(Exponent: LInt): Fp9Int;
    function IsZero:boolean;
    function IsOne:boolean;
    function IsASquare(Var Sqrt:Fp9int):boolean;
    function toPowerP:Fp9Int;
    function Sqr:Fp9Int;
    function MultiplyByGamma:Fp9Int;
    procedure SetToRandom;
    procedure SetToZero;
    procedure SetToOne;
    procedure SetTowerParams(TParam:PtrTowerParams18);
    procedure SetFromStrings(a,b,c:String);
  end;

  procedure _Add_FP9(const Left,Right:Fp9Int; var Result:Fp9Int);
  procedure _Sub_FP9(Const Left,Right:Fp9Int;var Result: Fp9Int);
  procedure _Mul_FP9(const Left,Right:Fp9Int;var Result:Fp9Int);
  procedure _Mul_FP3_FP9(const Left: Fp3Int; Right: Fp9Int; var Result: Fp9Int);
  procedure _Mul_FP_FP9(const Left: Lint; Right: Fp9Int; var Result: Fp9Int);
  procedure _Mul_FP9_By_V(const value: Fp9Int; var Result: Fp9Int);
  procedure _Sqr_FP9(const value: Fp9Int; var Result: Fp9Int);
  function _Equals_FP9(const Left,Right:Fp9Int):boolean;
  procedure _Inv_FP9(const value:Fp9Int;var Result:Fp9Int);
  procedure _Neg_FP9(const Value: Fp9Int; var Result: Fp9Int);
  procedure _Pow_FP9(const value:Fp9Int;Exponent:LInt;var Result:Fp9Int);
  function _FP9_To_DecimalString(const value:Fp9Int):string;
  function _FP9_To_HexString(const value:Fp9Int):string;
  procedure _FP9_From_Strings(const valueA,valueB,valueC:String;var Result:Fp9Int);
  procedure _Sparse2_Mul_FP9(const Left,Right:Fp9Int;var Result:Fp9Int);
  procedure _Sparse3_Mul_FP9(const Left,Right:Fp9Int;var Result:Fp9Int);
  procedure _Sparse1_Mul_FP9(const Left,Right:Fp9Int;var Result:Fp9Int);
  procedure _Sparse_Add_FP9(const Left,Right:Fp9Int; var Result:Fp9Int);
  procedure _Nomod_Sub_FP9(Const Left,Right:Fp9Int;var Result: Fp9Int);
  procedure _Pow_FP9_P(const value:Fp9Int;var Result:Fp9Int);
  procedure _Pow_FP9_P_i(const value:Fp9Int; power: integer;var Result:Fp9Int);
  function _Is_Fp9_ASquare(const value:Fp9Int;var Sqrt:Fp9Int):boolean;


implementation

{*******************************************************************************}
            ///      Procedures for FP9 Arithmetic
{*******************************************************************************}

    {**********   Add two Sparse FP9 integers *****************}
    //  Left is in Fp and Right.C is 0 (for D-type Twist only !)
procedure _Sparse_Add_FP9(const Left,Right:Fp9Int; var Result:Fp9Int);
begin
Result.Tower:=Left.Tower;
_Add_FP3(Left.a,Right.a,Result.a);
Result.b:=Right.b;
Result.c.a.Data.i16[-2]:=0;
Result.c.b.Data.i16[-2]:=0;
end;

    {**********   Add two FP9 integers *****************}
procedure _Add_FP9(const Left,Right:Fp9Int; var Result:Fp9Int);
begin
 Result.Tower:=Left.Tower;
_Add_FP3(Left.a,Right.a,Result.a);
_Add_FP3(Left.b,Right.b,Result.b);
_Add_FP3(Left.c,Right.c,Result.c);
end;

    {**********   Substract two FP9 integers ***********}
procedure _Sub_FP9(Const Left,Right:Fp9Int;var Result: Fp9Int);
begin
Result.Tower:=Left.Tower;
_Sub_FP3(Left.a,Right.a,Result.a);
_Sub_FP3(Left.b,Right.b,Result.b);
_Sub_FP3(Left.c,Right.c,Result.c);
end;

 {**********   Substract two FP9 integers ***********}
     // Without Modular reduction (for optimization in some cases)
procedure _Nomod_Sub_FP9(Const Left,Right:Fp9Int;var Result: Fp9Int);
begin
Result.Tower:=Left.Tower;
_Nomod_Sub_FP3(Left.a,Right.a,Result.a);
_Nomod_Sub_FP3(Left.b,Right.b,Result.b);
_Nomod_Sub_FP3(Left.c,Right.c,Result.c);
end;

    {********** Sparse Multiply two FP6 integers *****************}
    // Sparse_1 Multiplication Using Karatsuba
    // Right is in FP3 (Alwyas Sparse)  , for D-type Twist only !
procedure _Sparse1_Mul_FP9(const Left,Right:Fp9Int;var Result:Fp9Int);
begin
Result.Tower:=Left.Tower;
_Mul_FP3(Left.a,Right.a,Result.a);
_Mul_FP3(Left.b,Right.a,Result.b);
_Mul_FP3(Left.c,Right.a,Result.c);
end;

    {********** Sparse Multiply two FP6 integers *****************}
    //Sparse 1 Multiplication Using Karatsuba Multiplication
    // Right.C=0 and right.a=0 and right.b.b=0 (only right.b.a<>0). for M-type Twist only !
procedure _Sparse3_Mul_FP9(const Left,Right:Fp9Int;var Result:Fp9Int);
var tmp:array[0..1] of FP3Int;
begin
Result.Tower:=Left.Tower;
_Add_FP3(Left.b,Left.c,tmp[0]);
_Add_FP3(Left.a,Left.b,tmp[1]);
_Mul_FP3(Left.b,Right.b,Result.c);
_Mul_FP3(tmp[0],right.b,Result.a);
_Nomod_Sub_FP3(Result.a,Result.c,Result.a);
_Mul_FP3(Result.a,Left.Tower^.sigma,Result.a);
_Mul_FP3(tmp[1],right.b,Result.b);
_Sub_FP3(Result.b,Result.c,Result.b);
end;

        {********** Sparse Multiply two FP6 integers *****************}
        // Sparse-2 Multiplication Using Karatsuba Multiplication
        /// Right.c is always 0 (Partially Sparse)   , for D-type Twist only !)
procedure _Sparse2_Mul_FP9(const Left,Right:Fp9Int;var Result:Fp9Int);
var tmp:array[0..6] of FP3Int;
begin
Result.Tower:=Left.Tower;
_Mul_FP3(Left.a,Right.a,tmp[0]);
_Mul_FP3(Left.b,Right.b,tmp[1]);
_Nomod_Add_FP3(Left.b,Left.c,tmp[3]);
_NoMod_Add_FP3(Left.a,Left.b,tmp[4]);
_NoMod_Add_FP3(Right.a,Right.b,tmp[5]);
_Nomod_Add_FP3(Left.a,Left.c,tmp[6]);
_Mul_FP3(tmp[3],right.b,Result.a);
_Nomod_Sub_FP3(Result.a,tmp[1],Result.a);
_Mul_FP3(Result.a,Left.Tower^.sigma,Result.a);
_Add_FP3(tmp[0],Result.a,Result.a);
_Mul_FP3(tmp[4],tmp[5],Result.b);
_Nomod_Sub_FP3(Result.b,tmp[0],Result.b);
_Sub_FP3(Result.b,tmp[1],Result.b);
_Mul_FP3(tmp[6],right.a,Result.c);
_Nomod_Sub_FP3(Result.c,tmp[0],Result.c);
_Add_FP3(Result.c,tmp[1],Result.c);
end;



        {********** Multiply two FP9 integers *****************}
procedure _Mul_FP9(const Left,Right:Fp9Int;var Result:Fp9Int);
var tmp:array[0..9] of Fp3Int;
begin
        { Multiplication Using Karatsuba Multiplication  }
Result.Tower:=Left.Tower;
_Mul_FP3(Left.a,Right.a,tmp[0]);
_Mul_FP3(Left.b,Right.b,tmp[1]);
_Mul_FP3(Left.c,Right.c,tmp[2]);
_Add_FP3(Left.b,Left.c,tmp[3]);
_Add_FP3(Right.b,Right.c,tmp[4]);
_Add_FP3(Left.a,Left.b,tmp[5]);
_Add_FP3(Right.a,Right.b,tmp[6]);
_Add_FP3(Left.a,Left.c,tmp[7]);
_Add_FP3(Right.a,Right.c,tmp[8]);
_Mul_FP3(tmp[3],tmp[4],Result.a);
_Nomod_Sub_FP3(Result.a,tmp[1],Result.a);
_Sub_FP3(Result.a,tmp[2],Result.a);
_Mul_FP3(Result.a,Left.Tower^.sigma,Result.a);
_Add_FP3(tmp[0],Result.a,Result.a);
_Mul_FP3(tmp[5],tmp[6],Result.b);
_Nomod_Sub_FP3(Result.b,tmp[0],Result.b);
_Sub_FP3(Result.b,tmp[1],Result.b);
_Mul_FP3(tmp[2],Left.Tower^.sigma,tmp[9]);
_Add_FP3(Result.b,tmp[9],Result.b);
_Mul_FP3(tmp[7],tmp[8],Result.c);
_Nomod_Sub_FP3(Result.c,tmp[0],Result.c);
_Nomod_Sub_FP3(Result.c,tmp[2],Result.c);
_Add_FP3(Result.c,tmp[1],Result.c);
end;

        {********** Multiply FP3 with FP9 *****************}
procedure _Mul_FP3_FP9(const Left: Fp3Int; Right: Fp9Int; var Result: Fp9Int);
begin
Result.Tower:=Right.Tower;
_Mul_FP3(Right.a,Left,Result.a);
_Mul_FP3(Right.b,Left,Result.b);
_Mul_FP3(Right.c,Left,Result.c);
end;

        {********** Multiply FP with FP9 *****************}
procedure _Mul_FP_FP9(const Left: Lint; Right: Fp9Int; var Result: Fp9Int);
begin
Result.Tower:=Right.Tower;
_Mul_Fp_FP3(Left,Right.a,Result.a);
_Mul_FP_FP3(Left,Right.b,Result.b);
_Mul_FP_FP3(Left,Right.c,Result.c);
end;

        {********** Get Square of an FP9 *****************}
procedure _Sqr_FP9(const value: Fp9Int; var Result: Fp9Int);
var tmp: array [0..4] of Fp3Int;
begin
Result.Tower:=Value.Tower;
      {       Using Karatsuba Squering    }
_Sqr_FP3(Value.a,tmp[0]);
_Sqr_FP3(Value.b,tmp[1]);
_Sqr_FP3(Value.c,tmp[2]);
_Add_FP3(Value.a,Value.b,tmp[3]);
_Add_FP3(Value.a,Value.c,tmp[4]);
_Add_FP3(Value.b,Value.c,Result.a);
_Sqr_FP3(Result.a,Result.a);
_Sub_FP3(Result.a,tmp[1],Result.a);
_Sub_FP3(Result.a,tmp[2],Result.a);
_Mul_FP3(value.Tower^.Sigma,Result.a,Result.a);
_Add_FP3(Result.a,tmp[0],Result.a);
_Sqr_FP3(tmp[3],tmp[3]);
_Sub_FP3(tmp[3],tmp[0],tmp[3]);
_Sub_FP3(tmp[3],tmp[1],tmp[3]);
_Mul_FP3(Value.Tower^.Sigma,tmp[2],Result.b);
_Add_FP3(Result.b,tmp[3],Result.b);
_Sqr_FP3(tmp[4],Result.c);
_Sub_FP3(Result.c,tmp[0],Result.c);
_Sub_FP3(Result.c,tmp[2],Result.c);
_Add_FP3(Result.c,tmp[1],Result.c);
end;

        {**********   Compare two FP9 integers *************}
function _Equals_FP9(const Left,Right:Fp9Int):boolean;
begin
result:=(_Compare_FP3(Left.a,Right.a)=0) and (_Compare_FP3(Left.b,Right.b)=0)and (_Compare_FP3(Left.c,Right.c)=0);
end;

        {********** Get inverse of an FP9 *****************}
procedure _Inv_FP9(const value:Fp9Int;var Result:Fp9Int);
var t:array[0..9] of Fp3Int;
begin
Result.Tower:=Value.Tower;
_Sqr_FP3(value.a,t[0]);
_Sqr_FP3(value.b,t[1]);
_Sqr_FP3(value.c,t[2]);
_Mul_FP3(value.a,value.b,t[3]);
_Mul_FP3(value.a,value.c,t[4]);
_Mul_FP3(value.b,value.c,t[5]);
_Mul_FP3(value.Tower^.Sigma,t[5],t[7]);
_Sub_FP3(t[0],t[7],t[7]);
_Mul_FP3(value.Tower^.Sigma,t[2],t[8]);
_Sub_FP3(t[8],t[3],t[8]);
_Sub_FP3(t[1],t[4],t[9]);
_Mul_FP3(value.a,t[7],t[6]);
_Mul_FP3(value.Tower^.Sigma,value.c,t[0]);
_Mul_FP3(t[0],t[8],t[0]);
_Add_FP3(t[6],t[0],t[6]);
_Mul_FP3(value.Tower^.Sigma,value.b,t[0]);
_Mul_FP3(t[0],t[9],t[0]);
_Add_FP3(t[6],t[0],t[6]);
_Inv_FP3(t[6],t[6]);
_Mul_FP3(t[7],t[6],Result.a);
_Mul_FP3(t[8],t[6],Result.b);
_Mul_FP3(t[9],t[6],Result.c);
end;

       {****** Test if an Fp9 is a Sequare *******}
      // HAN, Dong-Guk, CHOI, Dooho, et KIM, Howon. Improved computation of square roots in specific finite fields. IEEE Transactions on Computers, 2009,
      //vol. 58, no 2, p. 188-196.
     // Adj, G., & Rodríguez-Henríquez, F. (2014). Square root computation over even extension fields. IEEE Transactions on Computers, 63(11), 2829-2841.
function _Is_Fp9_ASquare(const value:Fp9Int;var Sqrt:Fp9Int):boolean;
var a0,a1,a2:fp9Int;
    x,Beta:Lint;
begin
if value.IsZero then begin
                     result:=true;
                     Sqrt.SetToZero;
                     end
else begin
     if Value.Tower.FieldParam.p mod 4 =3 then begin  {p mod 4=3}
                                 _Sub_LInt(Value.Tower.FieldParam.p,3,x);
                                 _Shr_LInt(x,2);
                                 _pow_fp9(Value,x,a0);
                                 _Mul_FP9(a0,value,a1);
                                 _Sqr_FP9(a1,a1);
                                 _Mul_FP9(a1,a0,a1);
                                 _Pow_FP9_P(a0,a2);
                                 _Mul_FP9(a1,a2,a1);
                                 _Pow_FP9_P(a1,a1);

                                 _Pow_FP9_P(a1,a2);
                                 _Pow_FP9_P(a2,a2);
                                 _Mul_FP9(a1,a2,a1);
                                 _Pow_FP9_P(a2,a2);
                                 _Pow_FP9_P(a2,a2);
                                 _Mul_FP9(a1,a2,a1);
                                 _Pow_FP9_P(a2,a2);
                                 _Pow_FP9_P(a2,a2);
                                 _Mul_FP9(a1,a2,a1);

                                 _Mul_FP9(a1,a0,a1);
                                 _Mul_FP9(a1,value,Sqrt);

                                 _Mul_FP9(Sqrt,a1,a1);
                                 if a1.a.a=value.Tower.FieldParam.p-1 then Result:=false
                                 else result:=true;
                                 end
     else if Value.Tower.FieldParam.pmod8=5 then begin
                                                 _Sub_LInt(Value.Tower.FieldParam.p,5,x);
                                                 _Shr_LInt(x,3);
                                                 _Pow_Fp9(value,x,a0);
                                                 _Mul_FP9(a0,value,a1);

                                                 _Mul_FP9(a0,a1,a2);

                                                 _Sqr_FP9(a2,a2);
                                                 _Mul_FP9(a1,a2,a2);
                                                 _Pow_Fp9_P(a0,a1);
                                                 _Mul_FP9(a1,a2,a1);
                                                 _Pow_FP9_P(a1,a1);

                                                 _Pow_FP9_P(a1,a2);
                                                 _Pow_FP9_P(a2,a2);
                                                 _Mul_FP9(a1,a2,a1);
                                                 _Pow_FP9_P(a2,a2);
                                                 _Pow_FP9_P(a2,a2);
                                                 _Mul_FP9(a1,a2,a1);
                                                 _Pow_FP9_P(a2,a2);
                                                 _Pow_FP9_P(a2,a2);
                                                 _Mul_FP9(a1,a2,a1);
                                                 _Mul_FP9(a0,a1,a1);

                                                 _Sqr_FP9(a1,a2);
                                                 _Mul_FP9(a2,value,a2);
                                                  _Sqr_FP9(a2,a2);
                                                 if a2.a.a=value.Tower.FieldParam.p-1 then Result:=false
                                                 else result:=true;
                                                 _Mul_FP_Fp3(value.Tower.FieldParam.two_expq_div8, a1.a,a1.a);
                                                 _Mul_FP_Fp3(value.Tower.FieldParam.two_expq_div8, a1.b,a1.b);
                                                 _Mul_FP_Fp3(value.Tower.FieldParam.two_expq_div8, a1.c,a1.c);
                                                 _Mul_FP9(value,a1,a2);
                                                 _Mul_FP9(a2,a1,a0);
                                                 _Add_FP9(a0,a0,a0);
                                                 _Dec_LInt(a0.a.a,1);
                                                 _Mul_FP9(a2,a0,Sqrt);
                                                 end
     else result:=false;
     end;
end;

        {********** Get negative of an FP9 *****************}
procedure _Neg_FP9(const Value: Fp9Int; var Result: Fp9Int);
begin
Result.Tower:=Value.Tower;
_Neg_FP3(Value.a,Result.a);
_Neg_FP3(Value.b,Result.b);
_Neg_FP3(Value.c,Result.c);
end;

        {********** Multiply Fp9 by V *****************}
procedure _Mul_Fp9_By_V(const value: Fp9Int; var Result: Fp9Int);
var tmp:Fp3Int;
begin
if (Value.Tower.Sigma.a=0)and(Value.Tower.Sigma.c=0)and(Value.Tower.Sigma.b=1) then
begin
  tmp:=value.b;
  Result.Tower:=Value.Tower;
  Result.b:=value.a;
  _Mul_Fp3_u(Value.c,Result.a);
  Result.c:=tmp;
end
else begin
     Result.Tower:=Value.Tower;
     Result.b:=Value.a;
    Result.c:=Value.b;
    _Mul_FP3(Value.c,Value.Tower^.Sigma,Result.a);
     end;


end;

        {********** Raise an FP9 to a LInt power ************}
procedure _Pow_FP9(const value:Fp9Int;Exponent:LInt;var Result:Fp9Int);
var i:word;
    tmp:Fp9Int;
begin
Result.SetTowerParams(Value.Tower);
Result.a:=Value.a;
Result.b:=Value.b;
Result.c:=Value.c;
tmp:=Value;
if (Exponent>1) then
for i:=Exponent.BitLength-2 downto 0 do begin
                                        _Sqr_FP9(Result,Result);
                                        if _Is_BitSet_At(Exponent,i) then _Mul_FP9(Result,tmp,Result);
                                        end;
end;

        {********** Raise an Fp9 to " p " power ************}
        /// use Precomputed Frobenius Constants
procedure _Pow_FP9_P(const value:Fp9Int;var Result:Fp9Int);
begin
Result.SetTowerParams(value.Tower);
Result.a:=Value.a.toPowerP;
_Mul_FP3(Value.b.toPowerP,Value.Tower^.FrobeniusP_Const[0],Result.b);
_Mul_FP3(Value.c.toPowerP,Value.Tower^.FrobeniusP_Const[1],Result.c);
end;

        {********** Raise an Fp9 to " p^i " power ************}
        /// use Precomputed Frobenius Constant
procedure _Pow_FP9_P_i(const value:Fp9Int; power: integer;var Result:Fp9Int);
var i:integer;
begin
Result.SetTowerParams(Value.Tower);
Result := Value;
for i := 0 to power - 1 do _Pow_FP9_P(Result, Result);
end;

        {****** Convert an FP9 to a Decimal String **********}
function _FP9_To_DecimalString(const value:Fp9Int):string;
begin
if Value.c.IsZero then begin
                 if Value.b.IsZero then begin
                                  if Value.a.IsZero then Result:='0'
                                  else Result:=Value.a.ToDecimalString;
                                  end
                 else begin
                      if Value.b.IsOne then begin
                                      if Value.a.IsZero then Result:='v'
                                      else Result:='v +'+Value.a.ToDecimalString;
                                      end
                       else begin
                            if Value.a.IsZero then Result:='('+Value.b.ToDecimalString+') * v'
                            else Result:='('+Value.b.ToDecimalString+') * v +'+Value.a.ToDecimalString;
                            end;
                      end
                 end
else begin
     if Value.c.IsOne then begin
                     if Value.b.IsZero then begin
                                      if Value.a.IsZero then Result:='v^2'
                                      else Result:='v^2 + '+Value.a.ToDecimalString;
                                      end
                     else begin
                          if Value.b.IsOne then begin
                                          if Value.a.IsZero then Result:='v^2 + v'
                                          else Result:='v^2 + v +'+Value.a.ToDecimalString;
                                          end
                          else begin
                               if Value.a.IsZero then Result:='v^2 +('+Value.b.ToDecimalString+') * v'
                               else Result:='v^2 + ('+Value.b.ToDecimalString+') * v +'+Value.a.ToDecimalString;
                               end;
                         end
                    end
    else begin
         if Value.b.IsZero then begin
                          if Value.a.IsZero then Result:='('+Value.c.ToDecimalString+') * v^2'
                                      else Result:='('+Value.c.ToDecimalString+') * v^2 + '+Value.a.ToDecimalString;
                                      end
                          else begin
                               if Value.b.IsOne then begin
                                               if Value.a.IsZero then Result:='('+Value.c.ToDecimalString+') * v^2 +'+ 'v'
                                               else Result:='('+Value.c.ToDecimalString+') * v^2 +'+' v +'+Value.a.ToDecimalString;
                                               end
                               else begin
                                    if Value.a.IsZero then Result:='('+Value.c.ToDecimalString+') * v^2 +'+'('+Value.b.ToDecimalString+') * v'
                                    else Result:='('+Value.c.ToDecimalString+') * v^2 +'+'('+Value.b.ToDecimalString+') * v +'+Value.a.ToDecimalString;
                                    end;
                              end
                          end;
        end;
end;

        {****** Convert an FP9 to a Hexadecimal String *******}
function _FP9_To_HexString(const value:Fp9Int):string;
begin
if Value.c.IsZero then begin
                 if Value.b.IsZero then begin
                                  if Value.a.IsZero then Result:='0'
                                  else Result:=Value.a.ToHexString;
                                  end
                 else begin
                      if Value.b.IsOne then begin
                                      if Value.a.IsZero then Result:='v'
                                      else Result:='v +'+Value.a.ToHexString;
                                      end
                       else begin
                            if Value.a.IsZero then Result:='('+Value.b.ToHexString+') * v'
                            else Result:='('+Value.b.ToHexString+') * v +'+Value.a.ToHexString;
                            end;
                      end
                 end
else begin
     if Value.c.IsOne then begin
                     if Value.b.IsZero then begin
                                      if Value.a.IsZero then Result:='v^2'
                                      else Result:='v^2 + '+Value.a.ToHexString;
                                      end
                     else begin
                          if Value.b.IsOne then begin
                                          if Value.a.IsZero then Result:='v^2 + v'
                                          else Result:='v^2 + v +'+Value.a.ToHexString;
                                          end
                          else begin
                               if Value.a.IsZero then Result:='v^2 +('+Value.b.ToHexString+') * v'
                               else Result:='v^2 + ('+Value.b.ToHexString+') * v +'+Value.a.ToHexString;
                               end;
                         end
                    end
    else begin
         if Value.b.IsZero then begin
                          if Value.a.IsZero then Result:='('+Value.c.ToHexString+') * v^2'
                                      else Result:='('+Value.c.ToHexString+') * v^2 + '+Value.a.ToHexString;
                                      end
                          else begin
                               if Value.b.IsOne then begin
                                               if Value.a.IsZero then Result:='('+Value.c.ToHexString+') * v^2 +'+ 'v'
                                               else Result:='('+Value.c.ToHexString+') * v^2 +'+' v +'+Value.a.ToHexString;
                                               end
                               else begin
                                    if Value.a.IsZero then Result:='('+Value.c.ToHexString+') * v^2 +'+'('+Value.b.ToHexString+') * v'
                                    else Result:='('+Value.c.ToHexString+') * v^2 +'+'('+Value.b.ToHexString+') * v +'+Value.a.ToHexString;
                                    end;
                              end
                          end;
        end;
end;

        {****** Convert String to an FP9 (Decimal/Hex)*******}
procedure _FP9_From_Strings(const valueA,valueB,valueC:String;var Result:Fp9Int);
var i:integer;
    s1,s2,Val:string;
    valid:boolean;
begin
Result.a.SetFormString(valueA);
Result.b.SetFormString(valueB);
Result.c.SetFormString(valueC);
end;

{*******************************************************************************}
      ///  Definitions of an FP9 integer operators and functions
{*******************************************************************************}
class operator Fp9Int.Add(const Left, Right: Fp9Int): Fp9Int;
begin
_Add_FP9(Left,Right,Result);
end;

class operator Fp9Int.Subtract(const Left, Right: Fp9Int): Fp9Int;
begin
_Sub_FP9(Left,Right,Result);
end;

class operator Fp9Int.Multiply(const Left, Right: Fp9Int): Fp9Int;
begin
_Mul_FP9(Left,Right,Result);
end;

class operator Fp9Int.Multiply(const Left: Fp3Int;  Right: Fp9Int): Fp9Int;
begin
_Mul_FP3_FP9(Left,Right,Result);
end;

function Fp9Int.MultiplyByGamma: Fp9Int;
begin
_Mul_FP9_By_V(Self,Result);
end;

function Fp9Int.Sqr: Fp9Int;
begin
_Sqr_FP9(Self,Result);
end;

class operator Fp9Int.Equal(const Left, Right: Fp9Int): Boolean;
begin
result:=_Equals_FP9(Left,Right);
end;

function Fp9Int.Inverse: Fp9Int;
begin
_Inv_FP9(Self,Result);
end;

function Fp9Int.IsASquare(Var Sqrt:Fp9int): boolean;
begin
result:=_Is_Fp9_ASquare(Self,Sqrt)
end;

function Fp9Int.IsOne: boolean;
begin
Result:=(a.IsOne)and(b.IsZero)and(c.IsZero);
end;

function Fp9Int.IsZero: boolean;
begin
Result:=(a.IsZero)and(b.IsZero)and(c.IsZero);
end;

class operator Fp9Int.Negative(const Value: Fp9Int): Fp9Int;
begin
_Neg_FP9(Value,Result);
end;

class operator Fp9Int.NotEqual(const Left, Right: Fp9Int): Boolean;
begin
Result:=not _Equals_FP9(Left,Right)
end;

function Fp9Int.Pow(Exponent: LInt): Fp9Int;
begin
_Pow_FP9(Self,Exponent,Result);
end;

procedure Fp9Int.SetTowerParams(TParam:PtrTowerParams18);
begin
Tower:=TParam;
a.SetFieldParams(Tparam.FieldParam);
b.SetFieldParams(Tparam.FieldParam);
c.SetFieldParams(Tparam.FieldParam);
end;

procedure Fp9Int.SetToZero;
begin
a.SetToZero;
b.SetToZero;
c.SetToZero;
end;

procedure Fp9Int.SetFromStrings(a, b, c: String);
begin
_FP9_From_Strings(a,b,c,Self);
end;

procedure Fp9Int.SetToOne;
begin
a.SetToOne;
b.SetToZero;
c.SetToZero;
end;

procedure Fp9Int.SetToRandom;
begin
a.SetToRandom;
b.SetToRandom;
c.SetToRandom;
end;

function Fp9Int.toDecimalString: string;
begin
Result:=_FP9_To_DecimalString(Self)
end;

function Fp9Int.toHexString: string;
begin
Result:=_FP9_To_HexString(Self)
end;

function Fp9Int.toPowerP: Fp9Int;
begin
_Pow_FP9_P(Self,Result);
end;

end.
