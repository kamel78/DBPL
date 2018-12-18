unit LargeIntegers;


{********************************************************************************************
   Développer par Faraoun Kamel Mohamed
   Université Djilali Liabes -Sidi Bel Abbes - Algérie
   kamel_mh@yahoo.fr

  Fast arithmetic computation over multiprecision large integers, with corresponding computational
  algebric algorithmes (GCD,Inversion, Powering, Modular Powering, Primes generation,......

********************************************************************************************}

interface

uses System.Types,System.UITypes, System.SysUtils, System.Math,VCL.dialogs;
                //150
Const MaxLen=190 ;

  Const basediv2=2147483648;
      base=4294967296;
      base_1=4294967295;
      basearray:array [false..true] of Uint64=(base_1,base);
Type
       TInt64=packed record
         lo,hi:Dword;
       end;
       PLIntArrayForm=^LIntArrayForm;
       LIntArrayForm=array of ShortInt;
       PUint32=^UInt32;
        TData= record case integer of
           0:(i32:packed array [-1..MaxLen] of Dword);
           1:(i16:packed array [-2..MaxLen*2] of word);
           2:(i8:packed array [-4..MaxLen*4] of Byte);
           end;

        { Define a signed Big integer structure : -2^^(Maxlen*32)...+2^^(Maxlen*32)}

    PLint=^Lint;
    LInt=record
        Data:TData;
        id:char;
        Limit:Integer;
        InverseFordivision:PLInt;
        private
        public
        function ToDecimalString:string;
        function ToHexString:String;
        function ToNafArray:LIntArrayForm;
        function ToIntArray:LIntArrayForm;
        Function InversModulo(Modulo:LInt):LInt;
        Function Sqr:LInt;
        function BitLength:word;
        function TestBit(i:integer):boolean;
        function IsOdd:boolean;
        function IsEven:boolean;
        function IsZero:boolean;
        function Pow(e:longword):LInt;
        function PowMod(e, p: LInt;MgStr:Pointer=nil ): LInt;
        function IsASqrMod(p: LInt;MgStruct:Pointer): boolean;
        function Sqrt(p:LInt;MgStruct:Pointer): LInt;
        procedure SetToRandom(Limit:LInt;Seed:Word=0);
        procedure SetToRandomOnBits(nbits:Integer);
        function Absolute:LInt;
        class operator Add(const Left, Right: LInt): LInt;
        class operator Add(const Left:LInt; Right: word): LInt;
        class operator Add(const Left:Word; Right: LInt): LInt;
        class operator Subtract(const Left, Right: LInt): LInt;
        class operator Subtract(const Left: LInt;Right :Integer): LInt;
        class operator Multiply(const Left, Right: LInt): LInt;
        class operator Multiply(const Left: LInt; Right: Word): LInt;
        class operator Multiply(Left: integer; const Right: LInt): LInt;
        class operator Modulus(const Left: LInt; Right: LInt): LInt;
        class operator Divide(const Left: LInt; Right: LInt): LInt;
        class operator Implicit(const Value: string): LInt;
        class operator Implicit(const Value: Longint): LInt;
        class operator LeftShift(const Value: LInt; Shift: Integer): LInt;
        class operator RightShift(const Value: LInt; Shift: Integer): LInt;
        class operator Equal(const Left, Right: LInt): Boolean;
        class operator Equal(const Left:LInt; Right: DWord): Boolean;
        class operator NotEqual(const Left, Right: LInt): Boolean;
        class operator NotEqual(const Left:LInt; Right: Dword): Boolean;
        class operator GreaterThan(const Left, Right: LInt): Boolean;
        class operator GreaterThan(const Left:LInt; Right: Dword): Boolean;
        class operator GreaterThanOrEqual(const Left, Right: LInt): Boolean;
        class operator GreaterThanOrEqual(const Left:LInt; Right: Dword): Boolean;
        class operator LessThan(const Left, Right: LInt): Boolean;
        class operator LessThanOrEqual(const Left, Right: LInt): Boolean;
        class operator LessThan(const Left:LInt; Right: Dword): Boolean;
        class operator LessThanOrEqual(const Left:LInt; Right: Dword): Boolean;
        class operator Negative(const Value: LInt): LInt;
        end;

        LIntArray=Array of LInt;
        LintMatrix=array of LintArray;
        PLIntArray=^LIntArray;
        PtrMontgomeryStruct=^MontgomeryStruct;
        MontgomeryStruct=Record
                         modu,r,v,invr:LInt;
                         nbits:Dword;
                         end;


  PtrFieldParams=^FieldParams;
  FieldParams=record                      // Base Field parametres
    Beta:integer;                         // non-residus elements of the irredictible polynomial
    p:LInt;                               //Modulo of Fp
    p_1_div_2,pmod8,p_5_div_8,p_3_div_4,inv_2,sru,sru2,two_expq_div8,two_expq_1div27{For FP27 P^9 Forbenius}:LInt;                 // Precomputed values for Sqrt acceleration (p-1)/2 and 2^-1[mod] and Beta^((p-1)/6)
    MontgomeryData:PtrMontgomeryStruct;   // Montgomery parametres used for Exponontiation optimization
  end;

        procedure _Add_Lint(t1{eax},t2{edx}:Lint;var Result{ecx}:LInt);
        Procedure _Sub_LInt(t1{eax},t2{edx}:Lint;Var result{ecx}:LInt);
        Procedure _Mul_LInt(const t1,t2:Lint;var Result:LInt);
        Procedure _Div_Mod_LInt(a,b:LInt;var Q,res:LInt);
        Procedure _Mod_LInt(a,b:LInt;var res:LInt);

        procedure _Shl_LInt(var a:LInt;count:word);
        procedure _Shr_LInt(var a:LInt;const count:word);

        Procedure _Sqr_LInt(const a:Lint;Var result:LInt);
        Procedure _Pow_LInt(a:LInt;var  Result:LInt ;p:longword);
        Procedure _Mod_Pow_LInt(a,b,modulo:LInt;var Result:LInt);
        Procedure _Mg_Mod_Pow_LInt(a,b,n:LInt;MgStruct:PtrMontgomeryStruct;var Result:LInt);
        Procedure _Sqrt_Lint(Number:Lint;Var Result:Lint);

        Function _Compare_LInt(const t1,t2:LInt):Smallint;
        function _IsNull(a:LInt):boolean;
        function _IsOne(a:LInt):boolean;
        Function _IsNeg(a:LInt):boolean;
        Function _Str_To_LInt(s:string):LInt;
        Function _LInt_To_Str(a:LInt):String;
        Function _LInt_To_Hex(a:LInt):String;
        Function ArrayToHex(a:Tbytes):wideString;
        Function HexToArray(a:widestring):Tbytes;
        Function _Hex_To_LInt(a:String):LInt;
        Function ByteArrayToLint(a:TBytes):LInt;
        Function LIntToByteArray(a:LInt):Tbytes;
        Function _Num_Of_Bits_LInt(a:LInt):word;
        Function _Num_Of_Decimals_LInt(a:LInt):Word;
        Procedure _HCopy_LInt(a:Lint;var b:LInt);
        Function _Is_BitSet_At(a:LInt;i:integer):boolean;
        procedure _Set_LInt_BitAt(var a:LInt;i:integer;value:boolean);

        Procedure _Inc_LInt(Var a:LInt;p:Dword);
        Procedure _Dec_LInt(Var a:LInt;p:Dword);
        procedure _Inv_Mod_LInt(a,modulo:LInt;var Result:LInt);

        Procedure ExtendedEuclid(a,b:LInt;Var s,t:LInt);
        Function ResolveDiopantiene(a,b,c:LInt;var x,y:LInt):Boolean;
        Procedure GetRandomLIntOnBits(Var Result:LInt;Numbits:Integer);
        Procedure GetRandomLIntLowerThan(var Result:LInt;Limit:LInt;Seed:Word=0;SameBitlength:boolean=true);
        Procedure InitMontgomeryStruct(Modulo:LInt;var Struct:PtrMontgomeryStruct);
        Function IsLIntPrime(Number:LInt;Prob:Extended=0.9999):Boolean;
        function ModSquareLInt(a,modulo:LInt;Var Output:LInt;MgStruct:PtrMontgomeryStruct; test:boolean=true):boolean;
        procedure BinaryGCDLInt(a,b:LInt;var Result:LInt);
        function NAFToLInt(naf:LIntArrayForm):LInt;
        Function IsLIntOdd(a:LInt):boolean;
        procedure P_Add_Lint(const t1,t2:Lint;var result:Lint);
        procedure P_Sub_Lint(const t1,t2:Lint;var result:Lint);
        function HammingWeight(L:Lint;Negative:Boolean=False):Dword;
        function LIntToNAF(Val:LInt):LIntArrayForm;
        function LIntToIntArray(a:LInt):LIntArrayForm;
        function ArrayHammingWeight(L:LIntArrayForm):Dword;


 { $DEFINE PUREPASCAL}
             var numinvs:Longint;
implementation

{***************************************************************************************}
{             Addition of two Positive Large integeres                                  }
{***************************************************************************************}
procedure P_Add_Lint(const t1,t2:Lint;var result:Lint);
var i,LoopId,incarr,carr,val:DWord;
    BigIsFirst:boolean;
begin
BigIsFirst:=t1.Data.i16[-2]>t2.Data.i16[-2];
if  BigIsFirst then begin
                    LoopId:=t2.Data.i16[-2];
                    Result.Data.i16[-2]:=t1.Data.i16[-2];
                    end
else begin
     Loopid:=t1.Data.i16[-2];
     Result.Data.i16[-2]:=t2.Data.i16[-2];
     end;
Carr:=0;
if LoopId>0 then
for i:=0 to LoopId-1 do begin
                        val:=t1.Data.i32[i]+t2.Data.i32[i];
                        incarr:=Dword(val<t1.Data.i32[i]);
                        Result.Data.i32[i]:=val+carr;
                        carr:=incarr or Dword(Result.Data.i32[i]<carr);
                        end;
if BigIsFirst then begin
  if Result.Data.i16[-2]>0 then
  for i:=LoopId to Result.Data.i16[-2]-1 do begin
                                            Result.Data.i32[i]:=t1.Data.i32[i]+carr;
                                            carr:=Dword(Result.Data.i32[i] <carr);
                                            end
                   end
else
if Result.Data.i16[-2]>0 then
  for i:=LoopId to Result.Data.i16[-2]-1 do begin
                                            Result.Data.i32[i]:=t2.Data.i32[i]+carr;
                                            carr:=Dword(Result.Data.i32[i]<carr);
                                            end;
if Carr<>0 then begin
                inc(Result.Data.i16[-2]);
                Result.Data.i32[Result.Data.i16[-2]-1]:=Carr;
                end;
while(result.Data.i32[Result.Data.i16[-2]-1]=0)and(Result.Data.i16[-2]>0) do dec(Result.Data.i16[-2]);
end;
{***************************************************************************************}
{             Substraction of two Positive Large integeres                              }
{***************************************************************************************}
procedure P_Sub_Lint(const t1,t2:Lint;var result:Lint);
var i,LoopId,inborrow,borrow,val:DWord;
    BigIsFirst:boolean;
begin
BigIsFirst:=t1.Data.i16[-2]>=t2.Data.i16[-2];
if  BigIsFirst then begin
                    LoopId:=t2.Data.i16[-2];
                    Result.Data.i16[-2]:=t1.Data.i16[-2];
                    end
else begin
     Loopid:=t1.Data.i16[-2];
     Result.Data.i16[-2]:=t2.Data.i16[-2];
     end;
Borrow:=0;
if BigIsFirst then begin
                   if LoopId>0 then
                   for i:=0 to  LoopId-1 do begin
                                            val:=t1.Data.i32[i]-t2.Data.i32[i];
                                            inborrow:=Ord(val>t1.Data.i32[i]);
                                            Result.Data.i32[i]:=val-borrow;
                                            borrow:=inborrow or Ord(Result.Data.i32[i] =Dword(-1))and borrow;
                                            end;
                    if Result.Data.i16[-2]>0 then
                   for i:=LoopId to Result.Data.i16[-2]-1 do begin
                                                             Result.Data.i32[i]:=t1.Data.i32[i]-borrow;
                                                             borrow:=Ord(Result.Data.i32[i] =dword(-1)) and borrow;
                                                             end
                   end
else begin
     if Loopid>0 then
     for i:=0 to LoopId-1 do begin
                             val:=t2.Data.i32[i]-t1.Data.i32[i];
                             inborrow:=Ord(val>t2.Data.i32[i]);
                             Result.Data.i32[i]:=val-borrow;
                             borrow:=inborrow or Ord(Result.Data.i32[i] =Dword(-1))and borrow;
                             end;
     if Result.Data.i16[-2]>0 then
     for i:=LoopId to Result.Data.i16[-2]-1 do begin
                                               Result.Data.i32[i]:=t2.Data.i32[i]-Borrow;
                                               Borrow:=Ord(Result.Data.i32[i]= Dword(-1)) and borrow;
                                               end;

     end;
if borrow<>0 then begin
                  for i:=0 to Result.Data.i16[-2]-1 do begin
                                                       Result.Data.i32[i]:=not(Result.Data.i32[i]);
                                                       inborrow:=Ord(Result.Data.i32[i]=Dword(-1));
                                                       Result.Data.i32[i]:=Result.Data.i32[i]+borrow;
                                                       borrow:=Ord(Result.Data.i32[i]=Dword(-1))or inborrow;
                                                       end;
                  end;
while(result.Data.i32[Result.Data.i16[-2]-1]=0)and(Result.Data.i16[-2]>0) do dec(Result.Data.i16[-2]);
end;
{***************************************************************************************}
{             Addition of two signed Big integeres                                      }
{***************************************************************************************}
procedure _Add_Lint(t1{eax},t2{edx}:Lint;var result{ecx}:LInt);assembler;
{$IFDEF PUREPASCAL}
var tmp:word;
begin
if t1.Data.i16[-1]=t2.Data.i16[-1] then begin
                                        tmp:=t1.Data.i16[-1];
                                        P_Add_Lint(t1,t2,Result);
                                        Result.Data.i16[-1]:=tmp;
                                        end
else begin
     if _Compare_Lint(t1.Absolute,t2.Absolute)>=0 then tmp:=t1.Data.i16[-1]
     else tmp:=t2.Data.i16[-1];
     P_Sub_Lint(t1,t2,Result);
     Result.Data.i16[-1]:=tmp;
     end;
end;
{$ELSE}
{$IFDEF WIN32}
var carr,indic,zer:Dword;
    sizeA,sizeB:Word;
asm
         push ebx
         push esi
         push edi
         mov zer,0
         mov bx,Word ptr[edx]
         mov sizeB,bx
         mov bx,Word ptr [eax]   { Compare the size   }
         mov sizeA,bx
         shl sizeA,2
         shl sizeB,2
         cmp bx,Word ptr [edx]   { of the two numbers }
         jg @@1                  { and set the size   }
         mov bx,Word ptr [edx]   { of the result to   }
@@1:     mov Word ptr[ecx],bx    { the maximum one.   }
         cmp bx,0
         jne @033
         mov Word ptr[ecx+2],0
         jmp @@fin2
@033:    mov esi,ecx             { Save the result adress (ecx).}
         mov carr,0
         mov indic,0
         mov bx,word ptr [eax+2]{                        }
         cmp bx,Word ptr[edx+2] { Compare the signe of the two numbers}
         je @@same              { jump if same sign.                  }
        {--------------------------------------------------------}
         cmp bx,0                  {  put the negatif    }
         je @@t2                   {  number adr in edx  }
         xchg eax,edx              {                     }
         push ecx
         mov cx,sizeA
         mov bx,sizeB
         mov sizeA,bx
         mov sizeB,cx
         pop ecx
@@t2:    mov esi,ecx               { Save ecx.           }
         inc carr
         dec indic
        {--------------------------------------------------------}
@@same: mov word ptr[ecx+2],bx         {set the sign of the result (will be modified if diffrent signes))}
        xor ecx,ecx                   { Here start the addition with carry.  }
         mov cx,word ptr [esi]
         xor ebx,ebx
@@lp:    add ebx,4
         mov edi,0
         cmp bx,sizeB
         jg @@001
         mov edi,Dword ptr [edx+ebx]
@@001:   xor edi,indic                   {perform the 2 complement }
         add edi,carr
         mov carr,0
         adc carr,0
         cmp bx,sizeA
         jg @@002
         add edi,Dword ptr[eax+ebx]
         adc carr,0
@@002:   mov Dword ptr[esi+ebx],edi   { esi contain the initial ecx value (@t3)}
         or zer,edi
         loop @@lp
         mov ecx,carr
         cmp indic,0
         jne @@cmpRes             {Jump if the two numbers are with diffrents signs}
         add word ptr[esi],cx
         mov Dword ptr[esi+ebx+4],ecx
         or zer,ecx
         jmp @@fin
         {--------------------------------------------------------}
@@cmpRes:dec ecx
         mov word ptr[esi+2],cx
         cmp ecx,0
         je @@fin
         xor ecx,ecx
         mov cx,word ptr[esi]
         inc carr
         xor ebx,ebx
@@lp3:   add ebx,4
         mov edi,Dword ptr[esi+ebx]
         not edi
         add edi,carr
         mov Dword ptr[esi+ebx],edi
         mov carr,0
         adc carr,0
         loop @@lp3
         {--------------------------------------------------------}
@@fin:   cmp zer,0
         jne @@lll
         mov word ptr[esi],0
         jmp @@fin2
@@lll:   xor ecx,ecx
         mov cx,word ptr [esi]
         shl ecx,2
         add ecx,esi
@loop:   cmp dword ptr[ecx],0
         jne @@fin2
         dec Word ptr[esi]
         sub ecx,4
         jmp @loop
@@fin2:  mov ecx,esi
         pop edi
         pop esi
         pop ebx
end;
{$ELSE}
var carr,indic,zer:UInt64;
    sizeA,sizeB:Word;
asm
         push rbx
         mov zer,0
         xor rbx,rbx
         mov bx,Word ptr[RDX]
         shl bx,2
         mov dword ptr[rdx+rbx+4],0
         mov sizeB,bx
         mov bx,Word ptr[RDX]
         and bx,1
         jz @01
         add sizeB,4
@01:     mov bx,Word ptr [RCX]
         shl bx,2
         mov dword ptr[rcx+rbx+4],0
         mov sizeA,bx
         mov bx,Word ptr[RCX]
         and bx,1
         jz @02
         add sizeA,4
@02:     mov bx,Word ptr [RCX]
         cmp bx,Word ptr [RDX]
         jg @@1
         mov bx,Word ptr [RDX]
@@1:     mov Word ptr[R8],bx
		     cmp bx,0
         je @@resultnull
         jmp @@notnull
@@resultnull:mov Word ptr[R8+2],0
             jmp @@fin2
@@notnull:mov carr,0
         mov indic,0
         mov bx,word ptr [RCX+2]
		     cmp bx,Word ptr[RDX+2]
         je @@same
        {--------------------------------------------------------}
         cmp bx,0
         je @@t2
         xchg RCX,RDX
         mov ax,sizeA
         mov bx,sizeB
         mov sizeA,bx
         mov sizeB,ax
@@t2:    inc carr
         dec indic
        {--------------------------------------------------------}
@@same:  mov R9,RCX
         mov word ptr[R8+2],bx
         xor rcx,rcx
         mov cx,word ptr [R8]
         inc cx
         shr cx,1
         xor rbx,rbx
@@lp:    mov RAX,0
         cmp bx,sizeB
         jge @@001
         mov RAX,[RDX+RBX+4]
@@001:   xor RAX,indic
         add RAX,carr
         mov carr,0
         adc carr,0
         cmp bx,sizeA
         jge @@002
         add RAX,[R9+RBX+4]
         adc carr,0
@@002:   mov [R8+RBX+4],RAX
         or zer,RAX
         add RBX,8
         loop @@lp
         mov RCX,carr
         cmp indic,0
         jne @@cmpRes             {Jump if the two numbers are with diffrents signs}

         shr  rax,32
         mov edx,dword ptr[r8]
         and edx,1
         mul eax,edx
         test eax,eax
         setnz ch
         add cl,ch

         add byte ptr [R8],cl
         mov Qword ptr[R8+RBX+4],rcx
         or zer,RCX
         jmp @@fin
 {--------------------------------------------------------}
@@cmpRes:dec cx
         mov [R8+2],cx
         cmp cx,0
         je @@fin
         xor RCX,RCX
         mov cx,word ptr[R8]
         inc cx
         shr cx,1
         inc carr
         xor RBX,RBX
@@lp3:   mov RAX,[R8+RBX+4]
         not RAX
         add RAX,carr
         mov [R8+RBX+4],RAX
         mov carr,0
         adc carr,0
         add RBX,8
         loop @@lp3
        {--------------------------------------------------------}
@@fin:   cmp zer,0
         jne @@lll
         mov [R8],0
         jmp @@fin2
@@lll:   xor RCX,RCX
         mov cx,word ptr [R8]
         shl ecx,2
         add RCX,R8
@loop:   cmp dword ptr[RCX],0
         jne @@fin2
         dec Word ptr[R8]
         sub rcx,4
         jmp @loop
@@fin2:  pop rbx
end;
{$ENDIF}
{$ENDIF}

{***************************************************************************************}
{             Substraction of two signed Big integeres                                  }
{***************************************************************************************}
procedure _Sub_LInt(t1{eax},t2{edx}:Lint;var result{ecx}:LInt);assembler;
{$IFDEF PUREPASCAL}
var tmp:word;
begin
if t1.Data.i16[-1]<>t2.Data.i16[-1] then begin
                                         tmp:=t1.Data.i16[-1];
                                         P_Add_Lint(t1,t2,Result);
                                         Result.Data.i16[-1]:=tmp;
                                         end
else begin
     if _Compare_Lint(t1,t2)>=0 then tmp:=0
     else tmp:=$ffff;
     P_Sub_Lint(t1,t2,Result);
     Result.Data.i16[-1]:=tmp;
     end;
end;
{$ELSE}
{$IFDEF WIN32}
asm
        push esi
        mov esi,edx
        cmp eax,edx
        je @@zero             {compare if t1 and t2 are same }
        not word ptr[edx+2]   {invesre signe of t2}
        call _Add_Lint           {and perform an addition }
        cmp esi,ecx           {in the case of t2=Result don't inverse the result sign}
        je @@ret
        not word ptr[esi+2]
@@ret:  pop esi
        Ret
         {--------------------------------------------------------}
@@zero: push esi              { if t1 and t2 are same then return 0}
        mov esi,ecx
        mov ecx,DWord ptr[ecx]
        and ecx,$FFFFFFFF
        jz @@fin
@@lp:   mov Dword ptr[esi+ecx*4],0
        loop @@lp
        mov word ptr[esi],0
         {--------------------------------------------------------}
@@fin:  mov ecx,esi
        pop esi
end;
{$ELSE}
asm
        mov R10,RDX
        cmp RCX,RDX
        je @@zero             {compare if t1 and t2 are same }
        not word ptr[RDX+2]   {invesre signe of t2}
        call _Add_LInt           {and perform an addition }
        cmp R10,R8           {in the case of t2=Result don't inverse the result sign}
        je @@ret
        not word ptr[R10+2]
@@ret:  Ret
         {--------------------------------------------------------}
@@zero: mov [R8],0				{ if t1 and t2 are same then return 0}
@@fin:
end;
{$ENDIF}
{$ENDIF}


{***************************************************************************************}
{  Multiplication  of two signed Big integeres based on standard grade school algorithm }
{***************************************************************************************}
                {t1, t2 and Result Must be diffrent }
Procedure _Mul_LInt(const t1,t2:Lint;var Result:LInt);assembler;
{$IFDEF PUREPASCAL}
var i,j:integer;
    val:Uint64;
    Carr1,Carr2:Dword;
begin
if (t1.Data.i16[-2]=0)or(t2.Data.i16[-2]=0) then begin
                                                 Result.Data.i16[-2]:=0;
                                                 exit;
                                                 end;
for i:=0 to t2.Data.i16[-2]-1 do
            begin
            Carr1:=0;
            for j:=0 to t1.Data.i16[-2]-1 do begin
                                             if i=0 then Result.Data.i32[i+j]:=0;
                                             val:=Uint64(t2.Data.i32[i])*Uint64(t1.Data.i32[j])+Uint64(Carr1)+Uint64(Result.Data.i32[i+j]);
                                             Carr1:=Tint64(Val).hi;
                                             Result.Data.i32[i+j]:=Tint64(Val).lo;
                                             end;
            Result.Data.i32[i+t1.Data.i16[-2]]:=Carr1;
            end;
Result.Data.i16[-2]:=t1.Data.i16[-2]+t2.Data.i16[-2];
Result.Data.i16[-1]:=t1.Data.i16[-1] xor t2.Data.i16[-1];
while(result.Data.i32[Result.Data.i16[-2]-1]=0)and(Result.Data.i16[-2]>0) do dec(Result.Data.i16[-2]);
end;
{$ELSE}
{$IFDEF WIN64}
var carr:UInt64;
sizet1,sizet2:Dword;
asm
        push rbx
        push rsi
        push rdi
        mov r10,rdx
        mov bx,word ptr[rcx+2]
        xor bx,word ptr[r10+2]
        mov word ptr[r8+2],bx
        xor rbx,rbx
        xor rsi,rsi
        xor rdi,rdi
        mov esi,1
        mov bx,word ptr[rcx]
        mov sizet1,ebx
        {---------------------------------------}
        shl bx,2
        add rbx,4
        mov Dword ptr[RCX+RBX],0
        xor rbx,rbx
        {---------------------------------------}
        inc sizet1
        shr sizet1,1
        mov bx,word ptr[r10]
        mov sizet2,ebx
        {---------------------------------------}
        shl bx,2
        add rbx,4
        mov Dword ptr[r10+RBX],0
        xor rbx,rbx
        {---------------------------------------}
        inc sizet2
        shr sizet2,1
        cmp sizet2,0
        jz @isnull
        cmp sizet1,0
        jnz @start
@isnull:mov Dword ptr[r8],0             {one of the operand is null}
        jmp @@end4

@start: sub r8,4
        sub rcx,4
        sub r10,4
@@loop1:cmp esi,sizet1
        jg @@end1
        mov edi,1
        mov carr,0
@@loop2:cmp edi,sizet2
        jg @@end2
        xor rbx,rbx
        mov ebx,esi
        add ebx,edi
        dec ebx
        cmp esi,1
        jg @@next
        mov Qword ptr[r8+rbx*8],0
@@next: mov rax,Qword ptr[rcx+rsi*8]
        mul Qword ptr [r10+rdi*8]
        add rax,carr
        adc rdx,0
        add rax,Qword ptr[r8+rbx*8]
        adc rdx,0
        mov carr,rdx
        mov Qword ptr[r8+rbx*8],rax
        inc edi
        jmp @@loop2
@@end2: mov rbx,carr
        add edi,esi
        dec edi
        mov Qword ptr[r8+rdi*8],rbx
        inc esi
        jmp @@loop1
@@end1: mov bx,word ptr [rcx+4]
        add bx,word ptr [r10+4]
@@end3: mov word ptr[r8+4],bx
         {--------------------------------------------------------}
        add r8,4
        mov rsi,r8
        xor rcx,rcx
        mov cx,word ptr [rsi]
        shl ecx,2
        add rcx,rsi
@loop:  cmp dword ptr[rcx],0
        jne @@end4
        cmp rcx,rsi
        je @@end4
        dec Word ptr[rsi]
        sub rcx,4
        jmp @loop
         {--------------------------------------------------------}
@@end4: pop rdi
        pop rsi
        pop rbx
end;
{$ELSE}
var carr:Dword;
sizet1,sizet2:Dword;
asm
        push ebx
        push esi
        push edi
        mov bx,word ptr[eax+2]
        xor bx,word ptr[edx+2]
        mov word ptr[ecx+2],bx
        mov esi,1
        xor ebx,ebx
        mov bx,word ptr[eax]
        mov sizet1,ebx
        mov bx,word ptr[edx]
        mov sizet2,ebx
        cmp ebx,0
        jz @@isnull
        cmp sizet1,0
        jnz @@loop1
@@isnull:mov Dword ptr[ecx],0             {one of the operand is null}
        jmp @@end4
@@loop1:cmp esi,sizet1
        jg @@end1
        mov edi,1
        mov carr,0
@@loop2:cmp edi,sizet2
        jg @@end2
        mov ebx,esi
        add ebx,edi
        dec ebx
        cmp esi,1
        jg @@next
        mov Dword ptr[ecx+ebx*4],0
@@next: push eax
        push edx
        mov eax,Dword ptr[eax+esi*4]
        mul Dword ptr [edx+edi*4]
        add eax,carr
        adc edx,0
        add eax,Dword ptr[ecx+ebx*4]
        adc edx,0
        mov carr,edx
        mov Dword ptr[ecx+ebx*4],eax
        pop edx
        pop eax
        inc edi
        jmp @@loop2
@@end2: mov ebx,carr
        add edi,esi
        dec edi
        mov Dword ptr[ecx+edi*4],ebx
        inc esi
        jmp @@loop1
@@end1: mov bx,word ptr [eax]
        add bx,word ptr [edx]
@@end3: mov word ptr[ecx],bx
         {--------------------------------------------------------}
        mov esi,ecx
        xor ecx,ecx
        mov cx,word ptr [esi]
        shl ecx,2
        add ecx,esi
@loop:  cmp dword ptr[ecx],0
        jne @@end4
        cmp ecx,esi
        je @@end4
        dec Word ptr[esi]
        sub ecx,4
        jmp @loop
         {--------------------------------------------------------}
@@end4: pop edi
        pop esi
        pop ebx
end;
{$ENDIF}
{$ENDIF}

{***************************************************************************************}
{                       Compute the Square of a Big integer                             }
{***************************************************************************************}
Procedure _Sqr_LInt(const a:Lint;Var result:LInt);inline;
begin
_Mul_LInt(a,a,Result);
end;

{***************************************************************************************}
{             Fast division function based on Knuth division algorithm                  }
{***************************************************************************************}
Procedure _Div_Mod_LInt(a,b:LInt;var Q,res:LInt);inline;
{Const basediv2=2147483648;
      base=4294967296;
      base_1=4294967295;
      basearray:array [false..true] of Uint64=(base_1,base);}
var d,inv_d:Byte;
    bn,bn_1,ri,ri_1,ri_2:dword;
    i,j,hb,ha,hr,ptr:integer;
    rv,rhat,qhat,right,left,borrow,carry:Uint64;
    sigA,sigB:Word;
begin
if (b.Data.i16[-2]=0) then begin
                  Raise ERangeError.Create( 'Illegal Division by 0.');
                  exit;
                  end;
        {Determination of Q and Res signes}
q.data.i16[-1]:=a.data.i16[-1] and not(b.data.i16[-1]);
Res.data.i16[-1]:=a.data.i16[-1] or b.data.i16[-1];
        {a and b comparison  }
sigA:=a.data.i16[-1];
sigB:=b.data.i16[-1];   {save the signes of a and b}
a.data.i16[-1]:=0;
b.data.i16[-1]:=0;      {take the absolute value of a and b}
i:=_Compare_LInt(a,b);
a.data.i16[-1]:=sigA;
b.data.i16[-1]:=sigB;        {restore the signes of a and b}
if i<=0 then begin
             if i=0 then begin
                         res.data.i16[-2]:=0;
                         q.data.i16[-2]:=1;
                         q.data.i32[0]:=1;
                         exit;
                         end
             else begin
                  for i:=-1 to a.data.i16[-2]-1 do res.data.i32[i]:=a.data.i32[i];
                  q.data.i16[-2]:=0;
                  q.data.i32[0]:=0;
                  exit;
                  end;
             end;
ha:=a.data.i16[-2]-1;
//for i:=-1 to a.data.i16[-2]-1 do res.data.i32[i]:=a.data.i32[i];
Move(a.Data.i16[-2],res.Data.i16[-2],(a.data.i16[-2]+1)shl 2);
Res.data.i32[res.data.i16[-2]]:=0;
        { We can make a short division }
if (b.data.i16[-2]=1) then begin
                           rv:=0;
                           bn:=b.data.i32[0];
                           for i:=ha downto 0 do begin
                                                 rhat:=(rv shl 32)+res.data.i32[i];
                                                 Q.data.i32[i]:=rhat div bn;
                                                 rv:=rhat mod bn;
                                                 end;
                                     Q.data.i32[ha+1]:=0;
                                     Q.data.i16[-2]:=res.data.i16[-2];
                                     while (Q.data.i32[Q.data.i16[-2]-1]=0)and(Q.data.i16[-2]>0) do dec(Q.data.i16[-2]);
                                     if rv=0 then Res.data.i16[-2]:=0
                                     else begin
                                          Res.data.i16[-2]:=1;
                                          Res.data.i32[0]:=rv;
                                          end;
                                     exit;
                                     end;
        {   Structures initialization }
hb:=b.data.i16[-2]-1;
ptr:=ha-hb;
hr:=ha+1;
Q.data.i16[-2]:=ptr+3;
for i:=0 to Q.data.i16[-2] do Q.data.i32[i]:=0;
        {Camputation of the scaling factor d >[Base/Bn] }
bn:=b.data.i32[hb];
bn_1:=b.data.i32[hb-1];
d:=0;
while (bn<basediv2) do begin
                      bn:=bn shl 1;
                      inc(d);
                      end;
inv_d:=32-d;
        {Initialization of the divisor most signifients digits Bn and Bn-1 }
if d>0 then begin
            bn:=bn or(b.data.i32[hb-1] shr (inv_d));
            if hb>2 then bn_1:=(bn_1 shl d)or(b.data.i32[hb-2] shr(inv_d))
            else bn_1:=(bn_1 shl d);
            end;
                                 { The main Loop}
for i:=ptr  downto 0 do begin
                        {computation of the divident's most significant digits ri, ri_1 and ri_2 }
                        ri:=res.data.i32[hr];
                        ri_1:=res.data.i32[hr-1];
                        ri_2:=res.data.i32[hr-2];
                        if d>0 then begin
                                    ri:=(ri shl d)or(ri_1 shr(inv_d));
                                    ri_1:=(ri_1 shl d)or(ri_2 shr (inv_d));
                                    if (hr-3)<0 then ri_2:=ri_2 shl d
                                    else ri_2:=(ri_2 shl d)or(res.data.i32[hr-3] shr(inv_d));
                                    end;
                        {Estimation of the quotion digit qhat }
                        if ri<>bn then begin
                                       rhat:=(Uint64(ri) shl 32)or Uint64(ri_1);
                                       qhat:=rhat div bn;
                                       rhat:=rhat mod bn;
                                       right:=(rhat shl 32) or Uint64(ri_2);
                                       left:=bn_1*qhat;
                                       if left>right then begin
                                                          dec(qhat);
                                                          if (rhat+bn<base) and(left-bn_1>right+(Uint64(bn) shl 32)) then dec(qhat);
                                                          end;
                                       end
                        else begin
                             qhat:=base_1;
                             rhat:=Uint64(ri_1)+Uint64(bn);
                             right:=(rhat shl 32)or Uint64(ri_2);
                             if rhat<base then begin
                                               left:=bn_1*qhat;
                                               if left>right then begin
                                                                  dec(qhat);
                                                                  if (rhat+bn<base) and(left-bn_1>right+(Uint64(bn) shl 32)) then dec(qhat);
                                                                  end;
                                               end;
                             end;
                        {multiplication and substration  (aiai-1......ai-n):=(aiai-1......ai-n)-qhat*b  }
                        borrow:=base;
                        carry:=0;
                        for j:=i to i+hb do begin
                                            carry:=b.data.i32[j-i]*qhat+(carry shr 32);
                                            borrow:=Uint64(res.data.i32[j])+basearray[borrow>=base]-Uint64(carry and base_1);
                                            res.data.i32[j]:=borrow;
                                            end;
                        borrow:=Uint64(res.data.i32[i+hb+1])+basearray[borrow>=base]-Uint64((carry shr 32));
                        res.data.i32[i+hb+1]:=borrow;
                        {the quotion digit is stored }
                        Q.data.i32[i]:=qhat;
                        if (borrow<base) then begin
                                              carry:=0;
                                              for j:=i to i+hb do begin
                                                                  carry:=Uint64(res.data.i32[j])+Uint64(b.data.i32[j-i])+(carry shr 32);
                                                                  res.data.i32[j]:=carry;
                                                                  end;
                                              res.data.i32[i+hb+1]:=res.data.i32[j+1]+(carry shr 32);
                                              dec(Q.data.i32[i]);
                                              end;
                        hr:=hr-1;
                        end;
Res.data.i16[-2]:=b.data.i16[-2];
while (res.data.i32[res.data.i16[-2]-1]=0)and(Res.data.i16[-2]>0) do dec(Res.data.i16[-2]);
while (Q.data.i32[Q.data.i16[-2]-1]=0)and(Q.data.i16[-2]>0) do dec(Q.data.i16[-2]);
end;

{***************************************************************************************}
{             Fast modulo function based on Knuth division algorithm                    }
{***************************************************************************************}
Procedure _Mod_LInt(a,b:LInt;var res:LInt);
{$IFDEF WIN32}
Const basediv2=32768;
      base=65536;
      base_1=65535;
      basearray:array [false..true] of dword=(base_1,base);
var d,inv_d:Byte;
    bn,bn_1,ri,ri_1,ri_2:Word;
    i,j,hb,ha,hr,ptr:integer;
    rv,rhat,qhat,right,left,borrow,carry:Dword;
    sigA,sigB:Word;
    Negative:boolean;
    Q,Q2:Lint;
begin
if _IsNeg(a) and (a.Data.i16[-2]<>0) then begin
                 a.Data.i16[-1]:=0;
                 Negative:=true;
                 end
else Negative:=false;
if (b.Data.i16[-2]=0) then begin
                  Raise ERangeError.Create( 'Illegal Division by 0.');
                  exit;
                  end;
        {Determination of Res signes}
Res.data.i16[-1]:=a.data.i16[-1] or b.data.i16[-1];
        {a and b comparison  }
sigA:=a.data.i16[-1];
sigB:=b.data.i16[-1];   {save the signes of a and b}
a.data.i16[-1]:=0;
b.data.i16[-1]:=0;      {take the absolute value of a and b}
i:=_Compare_LInt(a,b);
a.data.i16[-1]:=sigA;
b.data.i16[-1]:=sigB;        {restore the signes of a and b}
if i<=0 then begin
             if i=0 then begin
                         res.data.i16[-2]:=0;
                         if Negative and (res.Data.i16[-2]<>0) then _Sub_LInt(b,res,Res);
                         exit;
                         end
             else begin
                  for i:=-1 to a.data.i16[-2]-1 do res.data.i32[i]:=a.data.i32[i];
                  if Negative and (res.Data.i16[-2]<>0) then _Sub_LInt(b,Res,Res);
                  exit;
                  end;
             end;

if (b.id='Y')and(a.Data.i16[-2]<=b.limit)then begin
                          _Mul_LInt(a,b.InverseFordivision^,Q);
                          _shr_lint(Q,b.limit shl 5);
                          _Mul_LInt(Q,b,Q2);
                          _Sub_LInt(a,Q2,Res);
                          if Negative  and (res.Data.i16[-2]<>0) then _Sub_LInt(b,res,Res);
                          if res.Data.i16[-1]<>0 then _add_LInt(b,res,Res);
                          exit;
                          end;
ha:=a.data.i16[-2]shl 1-1;
if a.data.i16[ha]=0 then dec(ha);
for i:=-1 to a.data.i16[-2]-1 do res.data.i32[i]:=a.data.i32[i];
Res.data.i32[res.data.i16[-2]]:=0;
        { We can make a short division }
if (b.data.i16[-2]=1)and(b.data.i16[1]=0) then begin
                                     rv:=0;
                                     bn:=b.data.i16[0];
                                     for i:=ha downto 0 do begin
                                                           rhat:=(rv shl 16)+res.data.i16[i];
                                                           rv:=rhat mod bn;
                                                           end;
                                     if rv=0 then Res.data.i16[-2]:=0
                                     else begin
                                          Res.data.i16[-2]:=1;
                                          Res.data.i32[0]:=rv;
                                          end;
                                     if Negative and not(res.Data.i16[-2]<>0) then _Sub_LInt(b,res,Res);
                                     exit;
                                     end;
        {   Structures initialization }
hb:=b.data.i16[-2]shl 1-1;
if b.data.i16[hb]=0 then dec(hb);
ptr:=ha-hb;
hr:=ha+1;
        {Camputation of the scaling factor d >[Base/Bn] }
bn:=b.data.i16[hb];
bn_1:=b.data.i16[hb-1];
d:=0;
while (bn<basediv2) do begin
                      bn:=bn shl 1;
                      inc(d);
                      end;
inv_d:=16-d;
        {Initialization of the divisor most signifients digits Bn and Bn-1 }
if d>0 then begin
            bn:=bn or(b.data.i16[hb-1] shr (inv_d));
            if hb>2 then bn_1:=(bn_1 shl d)or(b.data.i16[hb-2] shr(inv_d))
            else bn_1:=(bn_1 shl d);
            end;
                                 { The main Loop}
for i:=ptr  downto 0 do begin
                        {computation of the divident's most significant digits ri, ri_1 and ri_2 }
                        ri:=res.data.i16[hr];
                        ri_1:=res.data.i16[hr-1];
                        ri_2:=res.data.i16[hr-2];
                        if d>0 then begin
                                    ri:=(ri shl d)or(ri_1 shr(inv_d));
                                    ri_1:=(ri_1 shl d)or(ri_2 shr (inv_d));
                                    if (hr-3)<0 then ri_2:=ri_2 shl d
                                    else ri_2:=(ri_2 shl d)or(res.data.i16[hr-3] shr(inv_d));
                                    end;
                        {Estimation of the quotion digit qhat }
                        if ri<>bn then begin
                                       rhat:=(ri shl 16)or ri_1;
                                       qhat:=rhat div bn;
                                       rhat:=rhat mod bn;
                                       right:=(rhat shl 16) or ri_2;
                                       left:=bn_1*qhat;
                                       if left>right then begin
                                                          dec(qhat);
                                                          if (rhat+bn<base) and(left-bn_1>right+(bn shl 16)) then dec(qhat);
                                                          end;
                                       end
                        else begin
                             qhat:=base_1;
                             rhat:=ri_1+bn;
                             right:=(rhat shl 16)or ri_2;
                             if rhat<base then begin
                                               left:=bn_1*qhat;
                                               if left>right then begin
                                                                  dec(qhat);
                                                                  if (rhat+bn<base) and(left-bn_1>right+(bn shl 16)) then dec(qhat);
                                                                  end;
                                               end;
                             end;
                        {multiplication and substration  (aiai-1......ai-n):=(aiai-1......ai-n)-qhat*b  }
                        borrow:=base;
                        carry:=0;
                        for j:=i to i+hb do begin
                                            carry:=b.data.i16[j-i]*qhat+(carry shr 16);
                                            borrow:=res.data.i16[j]+basearray[borrow>=base]-(carry and base_1);
                                            res.data.i16[j]:=borrow;
                                            end;
                        borrow:=res.data.i16[i+hb+1]+basearray[borrow>=base]-(carry shr 16);
                        res.data.i16[i+hb+1]:=borrow;
                        if (borrow<base) then begin
                                              carry:=0;
                                              for j:=i to i+hb do begin
                                                                  carry:=res.data.i16[j]+b.data.i16[j-i]+(carry shr 16);
                                                                  res.data.i16[j]:=carry;
                                                                  end;
                                              res.data.i16[i+hb+1]:=res.data.i16[j+1]+(carry shr 16);
                                              end;
                        hr:=hr-1;
                        end;
Res.data.i16[-2]:=b.data.i16[-2];
if hb mod 2=0 then Res.data.i16[hb+1]:=0;
while (res.data.i32[res.data.i16[-2]-1]=0)and(Res.data.i16[-2]>0) do dec(Res.data.i16[-2]);
if Negative  and (res.Data.i16[-2]<>0) then _Sub_LInt(b,res,Res);

{$ELSE}

var d,inv_d:Byte;
    bn,bn_1,ri,ri_1,ri_2:dword;
    i,j,hb,ha,hr,ptr:integer;
    rv,rhat,qhat,right,left,borrow,carry:Uint64;
    sigA,sigB:Word;
	Negative:boolean;
  Q,Q2:Lint;
begin
if (a.Data.i16[-1]<>0) and (a.Data.i16[-2]<>0) then begin
                                                    a.Data.i16[-1]:=0;
                                                    Negative:=true;
                                                    end
else Negative:=false;
//if _isnull(b) then begin
if b.Data.i16[-2]=0 then begin
                         Raise ERangeError.Create( 'Illegal Division by 0.');
                         exit;
                         end;
        {Determination of Q and Res signes}
Res.data.i16[-1]:=a.data.i16[-1] or b.data.i16[-1];
        {a and b comparison  }
sigA:=a.data.i16[-1];
sigB:=b.data.i16[-1];   {save the signes of a and b}
a.data.i16[-1]:=0;
b.data.i16[-1]:=0;      {take the absolute value of a and b}
i:=_Compare_LInt(a,b);
a.data.i16[-1]:=sigA;
b.data.i16[-1]:=sigB;        {restore the signes of a and b}
if i<=0 then begin
             if i=0 then begin
                         res.data.i16[-2]:=0;
                         if Negative and (res.Data.i16[-2]<>0) then _Sub_LInt(b,res,Res);
                         exit;
                         end
             else begin
                  //for i:=-1 to a.data.i16[-2]-1 do res.data.i32[i]:=a.data.i32[i];
                  Move(a.Data.i16[-2],res.Data.i16[-2],(a.data.i16[-2]+1)shl 2);
                  if Negative and (res.Data.i16[-2]<>0) then _Sub_LInt(b,Res,Res);
                  exit;
                  end;
             end;
if (b.id='Y')and(a.Data.i16[-2]<=b.Limit)then begin
                          _Mul_LInt(a,b.InverseFordivision^,Q);
                          _shr_lint(Q,b.Limit shl 5);
                          //Move(Q.Data.i32[16],Q.Data.i32[0],64);
                          //Q.Data.i16[-2]:=Q.Data.i16[-2]-16;
                          //if Q.Data.i16[-2]<0 then Q.Data.i16[-2]:=0;
                          _Mul_LInt(Q,b,Q2);
                          _Sub_LInt(a,Q2,Res);
                          if res.Data.i16[-1]<>0 then _Add_Lint(b,res,res);
                          if Negative  and (res.Data.i16[-2]<>0) then _Sub_LInt(b,res,Res);
                          exit;
                          end;

ha:=a.data.i16[-2]-1;
//for i:=-1 to a.data.i16[-2]-1 do res.data.i32[i]:=a.data.i32[i];
Move(a.Data.i16[-2],res.Data.i16[-2],(a.data.i16[-2]+1)shl 2);
Res.data.i32[res.data.i16[-2]]:=0;
        { We can make a short division }
if (b.data.i16[-2]=1) then begin
                           rv:=0;
                           bn:=b.data.i32[0];
                           for i:=ha downto 0 do begin
                                                 rhat:=(rv shl 32)+res.data.i32[i];
                                                 rv:=rhat mod bn;
                                                 end;
                                     if rv=0 then Res.data.i16[-2]:=0
                                     else begin
                                          Res.data.i16[-2]:=1;
                                          Res.data.i32[0]:=rv;
                                          end;
									         if Negative and (res.Data.i16[-2]<>0) then _Sub_LInt(b,res,Res);
                                     exit;
                                     end;
        {   Structures initialization }
hb:=b.data.i16[-2]-1;

ptr:=ha-hb;
hr:=ha+1;
        {Camputation of the scaling factor d >[Base/Bn] }
bn:=b.data.i32[hb];
bn_1:=b.data.i32[hb-1];
d:=0;
while (bn<basediv2) do begin
                      bn:=bn shl 1;
                      inc(d);
                      end;
inv_d:=32-d;
        {Initialization of the divisor most signifients digits Bn and Bn-1 }
if d>0 then begin
            bn:=bn or(b.data.i32[hb-1] shr (inv_d));
            if hb>2 then bn_1:=(bn_1 shl d)or(b.data.i32[hb-2] shr(inv_d))
            else bn_1:=(bn_1 shl d);
            end;
                                 { The main Loop}
for i:=ptr  downto 0 do begin
                        {computation of the divident's most significant digits ri, ri_1 and ri_2 }
                        ri:=res.data.i32[hr];
                        ri_1:=res.data.i32[hr-1];
                        ri_2:=res.data.i32[hr-2];
                        if d>0 then begin
                                    ri:=(ri shl d)or(ri_1 shr(inv_d));
                                    ri_1:=(ri_1 shl d)or(ri_2 shr (inv_d));
                                    if (hr-3)<0 then ri_2:=ri_2 shl d
                                    else ri_2:=(ri_2 shl d)or(res.data.i32[hr-3] shr(inv_d));
                                    end;
                        {Estimation of the quotion digit qhat }
                        if ri<>bn then begin
                                       rhat:=(Uint64(ri) shl 32)or ri_1;
                                       qhat:=rhat div bn;
                                       rhat:=rhat mod bn;
                                       right:=(rhat shl 32) or ri_2;
                                       left:=bn_1*qhat;
                                       if left>right then begin
                                                          dec(qhat);
                                                          if (rhat+bn<base) and(left-bn_1>right+(Uint64(bn) shl 32)) then dec(qhat);
                                                          end;
                                       end
                        else begin
                             qhat:=base_1;
                             rhat:=Uint64(ri_1)+Uint64(bn);
                             right:=(rhat shl 32)or Uint64(ri_2);
                             if rhat<base then begin
                                               left:=bn_1*qhat;
                                               if left>right then begin
                                                                  dec(qhat);
                                                                  if (rhat+bn<base) and(left-bn_1>right+(Uint64(bn) shl 32)) then dec(qhat);
                                                                  end;
                                               end;
                             end;
                        {multiplication and substration  (aiai-1......ai-n):=(aiai-1......ai-n)-qhat*b  }
                        borrow:=base;
                        carry:=0;
                        for j:=i to i+hb do begin
                                            carry:=b.data.i32[j-i]*qhat+(carry shr 32);
                                            borrow:=Uint64(res.data.i32[j])+basearray[borrow>=base]-Uint64(carry and base_1);
                                            res.data.i32[j]:=borrow;
                                            end;
                        borrow:=res.data.i32[i+hb+1]+basearray[borrow>=base]-(carry shr 32);
                        res.data.i32[i+hb+1]:=borrow;
                        if (borrow<base) then begin
                                              carry:=0;
                                              for j:=i to i+hb do begin
                                                                  carry:=Uint64(res.data.i32[j])+Uint64(b.data.i32[j-i])+(carry shr 32);
                                                                  res.data.i32[j]:=carry;
                                                                  end;
                                              res.data.i32[i+hb+1]:=res.data.i32[j+1]+(carry shr 32);
                                              end;
                        hr:=hr-1;
                        end;
Res.data.i16[-2]:=b.data.i16[-2];
while (res.data.i32[res.data.i16[-2]-1]=0)and(Res.data.i16[-2]>0) do dec(Res.data.i16[-2]);
if Negative  and (res.Data.i16[-2]<>0) then _Sub_LInt(b,res,Res);
{$ENDIF}
end;

{***************************************************************************************}
{             Comparison of two signed Big integeres                                    }
{***************************************************************************************}
Function _Compare_LInt(const t1,t2:LInt):Smallint;
{Result         =  0:equal; 1:t1>t2; -1:t1<t2}
{$IFDEF PUREPASCAL}
var i:integer;
begin
if t1.Data.i16[-1]<>t2.Data.i16[-1] then Result:=t1.Data.i16[-1] or 1
else begin
     if t1.Data.i16[-2]>t2.Data.i16[-2] then result:=1 xor t1.Data.i16[-1]
     else if t1.Data.i16[-2]<t2.Data.i16[-2]  then result:=-1 xor t1.Data.i16[-1]
     else begin
          i:=t1.Data.i16[-2]-1;
          while (i>0)and(t1.Data.i32[i]=t2.Data.i32[i]) do i:=i-1;
          if t1.Data.i32[i]>t2.Data.i32[i] then result:=1 xor t1.Data.i16[-1]
          else if t1.Data.i32[i]<t2.Data.i32[i]  then result:=-1 xor t1.Data.i16[-1]
          else Result:=0;
          end;
     end;
end;
{$ELSE}
{$IFDEF WIN32}
asm
          push ebx
          push edi
          xor ebx,ebx
          xor cx,cx
          mov bx,word ptr[eax+2]
          cmp bx,Word ptr[edx+2]
          je @@samesig
          or bx,1
          mov cx,bx
          jmp @@end
@@samesig:xor cx,cx
          mov bx,word ptr[eax]
          cmp bx,word ptr[edx]
          je @@samsize
          sub bx,word ptr[edx]
          sets cl
          dec cx
          not cx
          xor cx,word ptr[eax+2]
          or cx,1
          jmp @@end
@@samsize:cmp bx,0
          je @@end
@@samsize1:mov edi,Dword ptr[eax+ebx*4]
          cmp edi,Dword ptr[edx+ebx*4]
          jne @@result
          dec bx
          jz @@end
          jmp @@samsize1
@@Result: setbe cl
          dec cx
          not cx
          or cx,1
          xor cx,word ptr[eax+2]
@@end:    mov result,cx
          pop edi
          pop ebx
end;
{$ELSE}
asm
          push rbx
          push rdi
          mov R9,RCX
          xor rbx,rbx
          xor rcx,rcx
          mov bx,word ptr[R9+2]
          cmp bx,Word ptr[RDX+2]
          je @@samesig
          or bx,1
          mov cx,bx
          jmp @@end
@@samesig:xor cx,cx
          mov bx,word ptr[R9]
          cmp bx,word ptr[RDX]
          je @@samsize
          sub bx,word ptr[RDX]
          sets cl
          dec cx
          not cx
          or cx,1
          xor cx,word ptr[R9+2]
          jmp @@end
@@samsize:cmp bx,0
          je @@end
@@samsize1:mov edi,Dword ptr[R9+rbx*4]
          cmp edi,Dword ptr[RDX+rbx*4]
          jne @@result
          dec bx
          jz @@end
          jmp @@samsize1
@@Result: setbe cl
          dec cx
          not cx
          or cx,1
          xor cx,word ptr[R9+2]
@@end:    mov result,cx
          pop rdi
          pop rbx
end;
{$ENDIF}
{$ENDIF}

{***************************************************************************************}
{             Right shift of a Big integer                                              }
{***************************************************************************************}
procedure _Shr_LInt(var a:LInt;const count:word);
{$IFDEF PUREPASCAL}
var i:integer;
    s,r:byte;
begin
s:=Count and 31;
r:=Count shr 5;
a.Data.i16[-2]:=a.Data.i16[-2]-r;
if r>0 then
    for i:=0 to a.Data.i16[-2] do a.Data.i32[i]:=a.Data.i32[i+r];
if s>0 then
    for i:=0 to a.Data.i16[-2]-2 do a.Data.i32[i]:=(a.Data.i32[i]shr s) or(a.Data.i32[i+1]shl (32-s));
a.Data.i32[a.Data.i16[-2]-1]:=a.Data.i32[a.Data.i16[-2]-1]shr s;
while(a.Data.i32[a.Data.i16[-2]-1]=0)and(a.Data.i16[-2]>0) do dec(a.Data.i16[-2]);
end;
{$ELSE}
{$IFDEF WIN32}
asm
       push esi
       push edi
       cmp Word ptr[eax],0
       je @@fin2
         {--------------------------------------------------------}
       xor ecx,ecx
       mov cx,Word ptr[eax]
       inc cx
       shl cx,2
       mov Dword ptr[eax+ecx],0
         {--------------------------------------------------------}
       xor ecx,ecx
       mov cx,count
       shr cx,5
       mov edi,eax
       jz @small
       add edi,4
       mov esi,edi
       //--------------- move bytes
       cmp cx,Word ptr[eax]
       jl @@lower
       mov cx,Word ptr[eax]
@@lower:sub Word ptr[eax],cx
       shl cx,2
       add esi,ecx
       mov cx,Word ptr[eax]
       rep movsd
       //----------------{rz with zeros}
       mov cx,Word ptr[eax]
       inc cx
       shl cx,2
       mov edi,eax
       add edi,ecx
       mov cx,count
       shr cx,5
       push eax
       xor eax,eax
       rep stosd
       pop eax
       mov edi,eax
       //---------------
@small:mov cx,count
       and cx,$1f
       jz @@fin2
       xor edx,edx
       mov dx,Word ptr[eax]
       shl edx,2
       add edx,eax
@@loop:cmp eax,edx
       jg @@fin
       add eax,4
       mov esi,Dword ptr[eax+4]
       shrd Dword ptr[eax],esi,cl
       mov esi,Dword ptr[eax]
       jmp @@loop
@@fin: cmp Dword ptr[eax-4],0
       jne @@fin2
       dec word ptr[edi]
@@fin2:pop edi
       pop esi
end;
{$ELSE}
asm
       push rsi
       push rdi
       cmp Word ptr[rcx],0
       je @@fin2
         {--------------------------------------------------------}
       mov r9,rcx
	   xor rcx,rcx
       mov cx,Word ptr[r9]
       inc cx
       shl cx,2
       mov Dword ptr[r9+rcx],0
         {--------------------------------------------------------}
       xor rcx,rcx
       mov cx,count
       shr cx,5
       mov rdi,r9
       jz @small
       add rdi,4
       mov rsi,rdi
       //--------------- move bytes
       cmp cx,Word ptr[r9]
       jl @@lower
       mov cx,Word ptr[r9]
@@lower:sub Word ptr[r9],cx
       shl cx,2
       add rsi,rcx
       mov cx,Word ptr[r9]
       rep movsd
       //----------------{rz with zeros}
       mov cx,Word ptr[r9]
       inc cx
       shl cx,2
       mov rdi,r9
       add rdi,rcx
       mov cx,count
       shr cx,5
       push rax
       xor rax,rax
       rep stosd
       pop rax
       mov rdi,r9
       //---------------
@small:mov cx,count
       and cx,$1f
       jz @@fin2
       xor rdx,rdx
       mov dx,Word ptr[r9]
       shl edx,2
       add rdx,r9
@@loop:cmp r9,rdx
       jg @@fin
       add r9,4
       mov esi,Dword ptr[r9+4]
       shrd Dword ptr[r9],esi,cl
       mov esi,Dword ptr[r9]
       jmp @@loop
@@fin: cmp Dword ptr[r9-4],0
       jne @@fin2
       dec word ptr[rdi]
@@fin2:pop rdi
       pop rsi
end;
{$ENDIF}
{$ENDIF}

{***************************************************************************************}
{              Left shift of a Big integer                                              }
{***************************************************************************************}
procedure _Shl_LInt(var a:LInt;count:word);
{$IFDEF PUREPASCAL}
var i:integer;
    s,r:byte;
begin
s:=Count and 31;
r:=Count shr 5;
a.Data.i16[-2]:=a.Data.i16[-2]+r;
if r>0 then
    for i:=a.Data.i16[-2] downto r do a.Data.i32[i]:=a.Data.i32[i-r];
    for i:=0 to r-1 do a.Data.i32[i]:=0;
if s>0 then begin
            if (a.Data.i32[a.Data.i16[-2]-1]shr (32-s))<>0 then begin
                                                           inc(a.Data.i16[-2]);
                                                           a.Data.i32[a.Data.i16[-2]-1]:=0;
                                                           end;
            for i:=a.Data.i16[-2]-1 downto 1 do a.Data.i32[i]:=(a.Data.i32[i]shl s) or(a.Data.i32[i-1]shr (32-s));
            a.Data.i32[0]:=a.Data.i32[0] shl s;
            end;
while(a.Data.i32[a.Data.i16[-2]-1]=0)and(a.Data.i16[-2]>0) do dec(a.Data.i16[-2]);
end;
{$ELSE}
{$IFDEF WIN32}
asm
        push esi
        push ebx
        push edi
        mov edi,eax
         {--------------------------------------------------------}
        xor ecx,ecx
        mov cx,Word ptr[eax]
        inc cx
        shl cx,2
        mov Dword ptr[eax+ecx],0
         {--------------------------------------------------------}
        xor ecx,ecx
        mov cx,count
        shr cx,5
        jz @@small
        add edi,4
        mov esi,edi
        //--------------- move bytes
        mov bx,cx
        mov cx,Word ptr[eax]
        shl cx,2
        add esi,ecx
        shr cx,2
        add cx,bx
        mov Word ptr[eax],cx
        shl cx,2
        add edi,ecx
        shr cx,2
        std
        rep movsd
        cld
       //----------------{rz with zeros}
        mov cx,count
        shr cx,5
        push eax
        xor eax,eax
        rep stosd
        pop eax
        mov edi,eax
        //---------------
@@small:mov cx,count
        and cx,$1f
        jz @@fin2
        xor edx,edx
        mov dx,Word ptr[eax]
       // and edx,maxlen-1
        shl edx,2
        add eax,edx
@@loop: cmp eax,edi
        je @@fin
        mov esi,Dword ptr[eax]
        shld Dword ptr[eax+4],esi,cl
        sub eax,4
        jmp @@loop
@@fin:  shl Dword ptr[edi+4],cl
        cmp Dword ptr[eax+edx+4],0
        je @@fin2
        inc Word ptr[edi]
@@fin2: pop edi
        pop ebx
        pop esi
end;
{$ELSE}
asm
       push rsi
        push rbx
        push rdi
        mov rdi,rcx
		mov r9,rcx
         {--------------------------------------------------------}
        xor rcx,rcx
        mov cx,Word ptr[rdi]
        inc cx
        shl cx,2
        mov Dword ptr[r9+rcx],0
         {--------------------------------------------------------}
        xor rcx,rcx
        mov cx,count
        shr cx,5
        jz @@small
        add rdi,4
        mov rsi,rdi
        //--------------- move bytes
        mov bx,cx
        mov cx,Word ptr[r9]
        shl cx,2
        add rsi,rcx
        shr cx,2
        add cx,bx
        mov Word ptr[r9],cx
        shl cx,2
        add rdi,rcx
        shr cx,2
        std
        rep movsd
        cld
       //----------------{rz with zeros}
        mov cx,count
        shr cx,5
        push rax
        xor eax,eax
        rep stosd
        pop rax
        mov rdi,r9
        //---------------
@@small:mov cx,count
        and cx,$1f
        jz @@fin2
        xor rdx,rdx
        mov dx,Word ptr[r9]
       // and edx,maxlen-1
        shl edx,2
        add r9,rdx
@@loop: cmp r9,rdi
        je @@fin
        mov esi,Dword ptr[r9]
        shld Dword ptr[r9+4],esi,cl
        sub r9,4
        jmp @@loop
@@fin:  shl Dword ptr[rdi+4],cl
        cmp Dword ptr[r9+rdx+4],0
        je @@fin2
        inc Word ptr[rdi]
@@fin2: pop rdi
        pop rbx
        pop rsi
end;
{$ENDIF}
{$ENDIF}

{***************************************************************************************}
{              Test if a Big integer is Null                                            }
{***************************************************************************************}
function _IsNull(a:LInt):boolean;inline;
begin
Result:=a.Data.i16[-2]=0;
end;

{***************************************************************************************}
{              Test if a Big integer is One                                            }
{***************************************************************************************}
function _IsOne(a:LInt):boolean;inline;
begin
Result:=(a.Data.i32[-1]=1)and(a.Data.i32[0]=1);
end;

{***************************************************************************************}
{              Test if a Big integer is Negative                                        }
{***************************************************************************************}
Function _IsNeg(a:LInt):boolean;
begin
result:=a.data.i16[-1]<>0;
end;

{***************************************************************************************}
{              Increment a Big Integer                                                  }
{***************************************************************************************}
Procedure _Inc_LInt(Var a:LInt;p:Dword);
{$IFDEF PUREPASCAL}
begin
_Add_Lint(a,p,a);
end;
{$ELSE}
{$IFDEF WIN32}
asm
        cmp word ptr[eax+2],0
        je @@start
        not word ptr[eax+2]
        call _Dec_LInt
        cmp word ptr[eax],0
        je @null
        not word ptr[eax+2]
@null:  Ret
@@start:push ebx
        xor ebx,ebx
        mov ecx,eax
        add ecx,4
        mov bx, word ptr[eax]
        cmp bx,0
        je @@endlp
@@non0: add Dword ptr[ecx],edx
        jnc @@end
        mov edx,1
        inc ebx
        shl ebx,2
        add ebx,eax
@@loop: add ecx,4
        cmp ecx,ebx
        je @@endlp
        add Dword ptr[ecx],1
        jnc @@end
        jmp @@loop
@@endlp:mov Dword ptr[ecx],edx
        inc word ptr[eax]
@@end:  pop ebx
end;
{$ELSE}
asm
        cmp word ptr[rcx+2],0
        je @@start
        not word ptr[rcx+2]
        call _Dec_LInt
        cmp word ptr[rcx],0
        je @null
        not word ptr[rcx+2]
@null:  Ret
@@start:push rbx
        xor rbx,rbx
        mov r9,rcx
        add r9,4
        mov bx, word ptr[rcx]
        cmp bx,0
        je @@endlp
@@non0: add Dword ptr[r9],edx
        jnc @@end
        mov edx,1
        inc ebx
        shl ebx,2
        add rbx,rcx
@@loop: add r9,4
        cmp r9,rbx
        je @@endlp
        add Dword ptr[r9],1
        jnc @@end
        jmp @@loop
@@endlp:mov Dword ptr[r9],edx
        inc word ptr[rcx]
@@end:  pop rbx
end;
{$ENDIF}
{$ENDIF}


{***************************************************************************************}
{              Decrement a Big Integer                                                  }
{***************************************************************************************}
Procedure _Dec_LInt(Var a:LInt;p:Dword);
{$IFDEF PUREPASCAL}
begin
_Sub_Lint(a,p,a);
end;
{$ELSE}
{$IFDEF WIN32}
asm
        cmp word ptr[eax+2],0
        je @@Start
        not word ptr[eax+2]
        call _Inc_LInt
        cmp word ptr[eax],0
        je @null
        not word ptr[eax+2]
@null:  Ret
@@start:push ebx
        push esi
        xor ebx,ebx
        mov ecx,eax
        add ecx,4
        mov bx, word ptr[eax]
        cmp bx,0
        je @@endlp
@@non0: mov esi,Dword ptr[ecx]
        sub esi,edx
        sub Dword ptr[ecx],edx
        jnb @@end
        mov edx,1
        inc ebx
        shl ebx,2
        add ebx,eax
@@loop: add ecx,4
        cmp ecx,ebx
        je @@endlp
        sub Dword ptr[ecx],1
        jnb @@end
        jmp @@loop
@@endlp:not Dword ptr[ecx-4]
        mov word ptr[eax+2],$ffff
@@end:  cmp Word ptr[eax],1
        jne @@fin1
        cmp DWord ptr[eax+4],0
        jne @@fin1
        mov Word ptr[eax],0
@@fin1: cmp Word ptr[eax],2
        jne @@fin
        cmp DWord ptr[eax+8],0
        jne @@fin
        mov Word ptr[eax],1
@@fin:  pop esi
        pop ebx
end;
{$ELSE}
asm
        cmp word ptr[rcx+2],0
        je @@Start
        not word ptr[rcx+2]
        call _Inc_LInt
        cmp word ptr[rcx],0
        je @Null
        not word ptr[rcx+2]
@null:  Ret
@@start:push rbx
        push rsi
        xor rbx,rbx
        mov r9,rcx
        add r9,4
        mov bx, word ptr[rcx]
        cmp bx,0
        je @@endlp
@@non0: mov esi,Dword ptr[r9]
        sub esi,edx
        sub Dword ptr[r9],edx
        jnb @@end
        mov edx,1
        inc ebx
        shl ebx,2
        add rbx,rcx
@@loop: add r9,4
        cmp r9,rbx
        je @@endlp
        sub Dword ptr[r9],1
        jnb @@end
        jmp @@loop
@@endlp:not Dword ptr[r9-4]
        mov word ptr[rcx+2],$ffff
@@end:  cmp Word ptr[rcx],1
        jne @@fin1
        cmp DWord ptr[rcx+4],0
        jne @@fin1
        mov Word ptr[rcx],0
@@fin1: cmp Word ptr[rcx],2
        jne @@fin
        cmp DWord ptr[rcx+8],0
        jne @@fin
        mov Word ptr[rcx],1
@@fin:  pop rsi
        pop rbx
end;
{$ENDIF}
{$ENDIF}


{***************************************************************************************}
{              Affect a Big integer to another                                          }
{***************************************************************************************}
Procedure _HCopy_LInt(a:Lint;var b:LInt);
{$IFDEF PUREPASCAL}
begin
Move(a.Data.i16[-2],b.Data.i16[-2],(a.data.i16[-2]+1)shl 2);
end;
{$ELSE}
{$IFDEF WIN32}
asm
    push esi
    push edi
    mov esi,eax
    mov edi,edx
    xor ecx,ecx
    mov cx,word ptr[eax]
    inc cx
    rep movsd
    pop edi
    pop esi
end;
{$ELSE}
asm
    push rsi
    push rdi
    mov rsi,rcx
    mov rdi,rdx
    xor ecx,ecx
    mov cx,word ptr[rsi]
    inc cx
    rep movsd
    pop rdi
    pop rsi
end;
{$ENDIF}
{$ENDIF}
{***************************************************************************************}
{             Find  the Integer Square root of a positive Large Integer                 }
{***************************************************************************************}
Procedure _Sqrt_Lint(Number:Lint;Var Result:Lint);
var a,a_plus1,b,tmp:Lint;
    nbit:Longint;
    stop:boolean;
begin
if (Number.Data.i16[-2]=0) then begin
                       Result.Data.i32[-1]:=0;
                       exit;
                       end
else if _IsNeg(Number) then begin
                           Raise ERangeError.Create( 'Can''t compute the square root of a negative number');
                           exit;
                           end;
nbit:=Number.BitLength shr 1;
a.data.i32[-1]:=nbit shr 5+1;
a.data.i32[a.data.i32[-1]-1]:=1 shl (nbit and 31);
Stop:=false;
while not Stop do begin
                  _Div_Mod_LInt(Number,a,b,tmp);
                  _Add_Lint(a,b,b);
                  _Shr_LInt(b,1);
                  _Sub_LInt(b,a,tmp);
                  tmp.data.i16[-1]:=0;
                  if _Compare_LInt(Tmp,Lint(2))=-1 then begin
                                                        _Mul_LInt(a,a,tmp);
                                                        if _Compare_LInt(tmp,Number)<1 then begin
                                                                                            _Add_Lint(a,1,a_plus1);
                                                                                            _Mul_LInt(a_plus1,a_plus1,tmp);
                                                                                            if _Compare_LInt(tmp,Number)=1 then begin
                                                                                                                                _HCopy_LInt(a,Result);
                                                                                                                                Stop:=true;
                                                                                                                                end;
                                                                                            end;
                                                        end;
                  _HCopy_LInt(b,a);
                  end;
end;
{***************************************************************************************}
{              Compute the Power of a Big integer to a longword value                   }
{***************************************************************************************}
Procedure _Pow_LInt(a:LInt;var  Result:LInt ;p:longword);
var i:byte;
    tmp:LInt;
begin
Result.data.i32[-1]:=1;
Result.data.i32[0]:=1;
for i:=31 downto 0 do begin
                      _Mul_LInt(result,result,tmp);
                      _HCopy_LInt(tmp,result);
                      if (p and(1 shl i)<>0)then begin
                                                 _Mul_LInt(Result,a,tmp);
                                                 _HCopy_LInt(tmp,Result);
                                                 end;
                      end;
end;

{***************************************************************************************}
{              Compute the Moular Power of a Big integer to a Big integer               }
{***************************************************************************************}
Procedure _Mod_Pow_LInt(a,b,modulo:LInt;var Result:LInt);
var i:word;
    tmp:LInt;
    nbits:Longint;
begin
Result.data.i32[-1]:=1;
Result.data.i32[0]:=1;
nbits:=_Num_Of_Bits_LInt(b);
for i:=nbits-1 downto 0 do begin
                         _Mul_LInt(result,result,tmp);
                         _Mod_LInt(tmp,modulo,Result);
                         if (b.data.i32[i shr 5] and (1 shl (i and 31))<>0) then begin
                                                                            _Mul_LInt(result,a,tmp);
                                                                            _Mod_LInt(tmp,modulo,Result);
                                                                            end;
                         end;
end;

{***************************************************************************************}
{              Test if the bit at position i is equal to 1  (i starts at 0)             }
{***************************************************************************************}
Function _Is_BitSet_At(a:LInt;i:integer):boolean;
begin
if ((i shr 5)>a.data.i16[-2]-1)  then result :=false
else  result:=a.data.i32[i shr 5] and (1 shl (i and 31))<>0;
end;

{***************************************************************************************}
{              Set a bit at position i to 1 or to 0                                     }
{***************************************************************************************}
procedure _Set_LInt_BitAt(var a:LInt;i:integer;value:boolean);
begin
if value then a.data.i32[i shr 5]:=a.data.i32[i shr 5] or (1 shl (i and 31))
else a.data.i32[i shr 5]:=a.data.i32[i shr 5] and (not(1 shl (i and 31)));
end;

{***************************************************************************************}
{              Return the number of bits of a Big integer                               }
{***************************************************************************************}
Function _Num_Of_Bits_LInt(a:LInt):word;
var l:longint;
    d:Dword;
begin
if a.data.i16[-2]=0 then Result:=0
else begin
     l:=(a.data.i16[-2]-1)shl 5;
     D:=a.data.i32[a.data.i16[-2]-1];
     While (D>0) do begin
                    inc(l);
                    d:=d shr 1;
                    end;
     Result:=l;
     end;
end;

{***************************************************************************************}
{              Return the number of decimal digits  of a Big integer                    }
{***************************************************************************************}
Function _Num_Of_Decimals_LInt(a:LInt):Word;
begin
Result:=Length(_LInt_to_Str(a));
end;

{***************************************************************************************}
{    Applay the exetended Euclid Algorithm to find s and t such That as+bt=gcd(a,b)     }
{***************************************************************************************}
Procedure ExtendedEuclid(a,b:LInt;Var s,t:LInt);inline;
var s0,s1,t0,t1,r,q,aa,bb,tmp:LInt;
begin
_HCopy_LInt(a,aa);
_HCopy_LInt(b,bb);
S0:=1;
s1.data.i32[-1]:=0;
t0.data.i32[-1]:=0;
t1:=1;
_Div_Mod_LInt(aa,bb,q,r);
While r.Data.i16[-2]<>0 do begin
                       _Mul_LInt(q,s1,tmp);
                       _Sub_LInt(s0,tmp,s);
                       _Mul_LInt(q,t1,tmp);
                       _Sub_LInt(t0,tmp,t);
                       _HCopy_LInt(s1,s0);
                       _HCopy_LInt(s,s1);
                       _HCopy_LInt(t1,t0);
                       _HCopy_LInt(t,t1);
                       _HCopy_LInt(bb,aa);
                       _HCopy_LInt(r,bb);
                       _Div_Mod_LInt(aa,bb,q,r);
                       end;
end;

{***************************************************************************************}
{             Resolve the linear Diopantiene Equation ax+by=c                           }
{***************************************************************************************}
Function ResolveDiopantiene(a,b,c:LInt;var x,y:LInt):Boolean;
var g,r,q,tmp:LInt;
begin
BinaryGCDLInt(a,b,g);
_Div_Mod_LInt(c,g,q,r);
if r.Data.i16[-2]=0 then begin
                  ExtendedEuclid(a,b,x,y);
                  _Mul_LInt(x,q,tmp);
                  _HCopy_LInt(tmp,x);
                  _Mul_LInt(y,q,tmp);
                  _HCopy_LInt(tmp,y);
                  Result:=True;
                  end
else begin
     Result:=False;
     Raise ERangeError.Create( 'No integer solution for this equation   ')
     end;
end;

{***************************************************************************************}
{                Compute the modular inverse of a Big integer                           }
{***************************************************************************************}
procedure _Inv_Mod_LInt(a,modulo:LInt;var Result:LInt);inline;  /// the modulo is supposed to be a prime
var Q,R,tmp,_a,_mod,s0,s1:LInt;
begin
inc(numinvs);
if a.Data.i16[-1]<>0 then _Mod_LInt(a,modulo,a);
_HCopy_LInt(Modulo,_mod);
_HCopy_LInt(a,r);
s0.Data.i32[-1]:=1;s0.Data.i32[0]:=1;
s1.data.i32[-1]:=0;
Q.Data.i32[-1]:=0;
While r.Data.i16[-2]<>0 do begin
                        _Mul_LInt(Q,s1,tmp);
                        _Sub_LInt(s0,tmp,Result);
                         _HCopy_LInt(s1,s0);
                         _HCopy_LInt(Result,S1);
                        _HCopy_LInt(_mod,_a);
                        _HCopy_LInt(R,_mod);
                        _Div_Mod_LInt(_a,_mod,Q,R);
                        end;
if Result.Data.i16[-1]<>0 then _Add_LInt(Result,modulo,Result);
end;


{***************************************************************************************}
{             Generate a random binary Big integer on Numbit Bits                       }
{***************************************************************************************}
Procedure GetRandomLIntOnBits(Var Result:LInt;Numbits:Integer);
var b,i:integer;
    x:Dword;
begin;
b:=Numbits shr 5;
Result.data.i32[-1]:=b;
Randomize;
For i:=0 to b-1 do Result.data.i32[i]:=Random(Maxlongint);
if Numbits and 31=0 then Result.data.i32[b-1]:=Result.data.i32[b-1] or $80000000
else begin
     inc(Result.data.i32[-1]);
     x:=(1 shl (Numbits and 31)-1);
     Result.data.i32[b]:=Random(x)or (1 shl ((Numbits and 31)-1));
     end;
end;

{***************************************************************************************}
{             Generate a random binary block Lower than Limit                           }
{***************************************************************************************}
Procedure GetRandomLIntLowerThan(var Result:LInt;Limit:LInt;Seed:Word=0;SameBitlength:boolean=true);
var i:integer;
begin
repeat
      if seed=0 then Randomize
      else RandSeed:=seed;
      Result.data.i32[-1]:=Limit.data.i32[-1];
      For i:=0 to Result.data.i32[-1]-1 do Result.data.i32[i]:=Random(MaxLongint);
      _Mod_LInt(Result,Limit,Result);
      while (Result.data.i32[Result.data.i16[-2]-1]=0)and(Result.data.i16[-2]>0) do dec(Result.data.i16[-2]);
      until (not(SameBitlength)or(Result.BitLength=Limit.BitLength));
end;

{***************************************************************************************}
{       Initialization of a montgomery strcutur with repect to a modulo (must be odd)   }
{***************************************************************************************}
Procedure InitMontgomeryStruct(Modulo:LInt;var Struct:PtrMontgomeryStruct);
var tmp,Head:DWord;
    i:integer;
begin
if (Modulo.data.i32[0] and 1=0) then begin
                                Raise ERangeError.Create( 'Montgomery exponentiation can operate on odd modulo only ');
                                exit;
                                end;
with Struct^ do begin
                 r.data.i32[-1]:=Modulo.data.i32[-1];
                 for i:=0 to r.data.i16[-2]-1 do r.data.i32[i]:=0;
                 Head:=Modulo.data.i32[Modulo.data.i16[-2]-1];
                 if Head>=$80000000 then begin
                                        nbits:=r.data.i16[-2]shl 5;
                                        inc(r.data.i16[-2]);
                                        r.data.i32[r.data.i16[-2]-1]:=1;
                                        end
                 else begin
                      tmp:=1;
                      nbits:=0;
                      while (tmp<=Head) do begin
                                           tmp:=tmp shl 1;
                                           inc(nbits);
                                           end;
                      r.data.i32[r.data.i16[-2]-1]:=tmp;
                      nbits:=nbits+(r.data.i16[-2]-1) shl 5;
                      end;
                 ResolveDiopantiene(r,modulo,1,Invr,v);
                 if _IsNeg(Invr) then _Add_LInt(Invr,Modulo,Invr);
                 v.data.i16[-1]:=not(v.data.i16[-1]);
                 if _IsNeg(v) then _Add_LInt(v,r,v);
                 _HCopy_LInt(Modulo,Modu);
                 end;
end;

{***************************************************************************************}
{          Compute a^^b (mod n) using montgomery algorithm (n must be odd number)       }
{               the montgomery approache is faster for very big numbers                 }
{***************************************************************************************}
Procedure _Mg_Mod_Pow_LInt(a,b,n:LInt;MgStruct:PtrMontgomeryStruct;var Result:LInt);
var s,r0,m,t,fia:LInt;
    i:integer;
    nbits:Longint;
begin
_HCopy_LInt(a,fia);
_Shl_LInt(fia,MgStruct^.nbits);
_Mod_LInt(fia,n,fia);
_HCopy_LInt(MgStruct^.r,r0);
nbits:=_Num_Of_Bits_LInt(b);
For i:=nbits-1 downto 0 do begin
                         _Mul_LInt(r0,r0,s);
                         _Mul_LInt(s,MgStruct^.v,t);
                         if t.data.i16[-2]>=MgStruct^.r.data.i16[-2] then begin
                                                               t.data.i16[-2]:=MgStruct^.r.data.i16[-2];
                                                               t.data.i32[t.data.i16[-2]-1]:=t.data.i32[t.data.i16[-2]-1] and (MgStruct^.r.data.i32[MgStruct^.r.data.i16[-2]-1]-1);
                                                               end;
                         while (t.data.i32[t.data.i16[-2]-1]=0)and(t.data.i16[-2]>0)do dec(t.data.i16[-2]);
                         _Mul_LInt(t,n,m);
                         _Add_LInt(s,m,r0);
                         _Shr_LInt(r0,MgStruct^.nbits);
                         if _Compare_LInt(r0,n)<>-1 then _Sub_LInt(r0,n,r0);
                         if (b.data.i32[i shr 5])and(1 shl (i and 31))<>0 then begin
                                                                          _Mul_LInt(r0,fia,s);
                                                                          _Mul_LInt(s,MgStruct^.v,t);
                                                                          if t.data.i16[-2]>=MgStruct^.r.data.i16[-2] then begin
                                                                                                                t.data.i16[-2]:=MgStruct^.r.data.i16[-2];
                                                                                                                t.data.i32[t.data.i16[-2]-1]:=t.data.i32[t.data.i16[-2]-1] and (MgStruct.r.data.i32[MgStruct.r.data.i16[-2]-1]-1);
                                                                                                                end;
                                                                          while (t.data.i32[t.data.i16[-2]-1]=0)and(t.data.i16[-2]>0)do dec(t.data.i16[-2]);
                                                                          _Mul_LInt(t,n,m);
                                                                          _Add_LInt(s,m,r0);
                                                                          _Shr_LInt(r0,MgStruct^.nbits);
                                                                          if _Compare_LInt(r0,n)<>-1 then _Sub_LInt(r0,n,r0);
                                                                          end;
                         end;
_Mul_LInt(r0,MgStruct^.invr,Result);
_Mod_LInt(Result,n,Result);
end;

{***************************************************************************************}
{       Test If a Big Integer is a prime number Using the Rabin-Miller Test             }
{***************************************************************************************}
Function IsLIntPrime(Number:LInt;Prob:Extended=0.9999):Boolean;
Var SminusOne,Base,x,Q,minusOne:LInt;
    nbItr,Rest:Byte;
    i:integer;
    nbits:Longint;
    Prim:Boolean;
    MgStr:PtrMontgomeryStruct;
begin
Randomize;
if (Number.data.i16[-2]=1)and(Number.data.i32[0]=2) then begin
                                                {The number is 2}
                                                Result:=True;
                                                exit;
                                                end;
if (Number.data.i32[0] and 1=0)  then begin
                                 {The number is Even}
                                 Result:=False;
                                 exit;
                                 end;
_Sub_LInt(Number,1,minusOne);

_Mod_Pow_LInt(2,MinusOne,Number,x);
if _Compare_LInt(x,1)<>0 then begin
                             Result:=false;
                             exit;
                             end;

_HCopy_LInt(MinusOne,SMinusOne);
_Shr_LInt(SminusOne,1);
nbItr:=Round((-Ln(1-Prob)/(2*ln(2))));
i:=0;
Prim:=True;
nbits:=_Num_Of_Bits_LInt(Number);
new(MgStr);
InitMontgomeryStruct(Number,MgStr);
while (i<=nbItr-1)and(Prim) do begin
                             GetRandomLIntOnBits(Base,nbits-1);
//                             base:='9519945929610798397375028530945473253896927751573450403237889582777';
                             _Mg_Mod_Pow_LInt(Base,minusOne,Number,MgStr,x);
                             if (x.data.i16[-2]<>1)or(x.data.i32[0]<>1) then Prim:=false
                             else begin
                                  _HCopy_LInt(SMinusOne,Q);
                                  Rest:=0;
                                  While (Rest=0)and Prim do begin
                                                            _Mg_Mod_Pow_LInt(Base,Q,Number,MgStr,x);
                                                            if _Compare_LInt(x,minusOne)=0 then Rest:=1
                                                            else begin
                                                                 if (x.data.i16[-2]<>1)or(x.data.i32[0]<>1) then
                                                                 Prim:=False
                                                                 else begin
                                                                      Rest:=Q.data.i32[0] and 1;
                                                                      _Shr_LInt(Q,1);
                                                                      end;
                                                                 end;
                                                            end;
                                  end;
                             i:=i+1;
                             end;
Result:=Prim;
end;

{***************************************************************************************}
{    Find  the Square of  a Big integer "a" modulo "modulo"  (Result^2==a[Modulo])      }
{***************************************************************************************}
function ModSquareLInt(a,modulo:LInt;Var Output:LInt;MgStruct:PtrMontgomeryStruct; test:boolean=true):boolean;inline;
var e1,e2,g,MinusOne,MinusOneDiv2,tmp,tmp2,tmp3:LInt;
begin
Result:=true;
if test then begin
                        {First we test if Modu is a Prime numbre}
             if (not IsLIntPrime(modulo)) then begin
                                              Result:=false;   //Raise ERangeError.Create( 'The modulo must be a prime number ');
                                              exit;
                                              end;
             end;
        {Then we test if a is a quadratic Residu modulo modu using Euler's Criterion}
_Sub_LInt(Modulo,1,tmp);
_Shr_LInt(tmp,1);
_Mg_Mod_Pow_LInt(a,tmp,modulo,MgStruct,g);
if (_Compare_LInt(a,0)<>0)and(_Compare_LInt(g,1)<>0) then begin
                                                    Result:=false; //Raise ERangeError.Create( _LInt_To_Str(a)+' is not a Quadratic residu Modulo '+_LInt_To_Str(modulo));
                                                    exit;
                                                    end;

if (Modulo.data.i32[0] and 3)=3 then begin
                                        {The modulo is in the form 4k+3 so direct computation...}
                               _HCopy_LInt(Modulo,Tmp);
                               _inc_lint(tmp,1);
                               _Shr_LInt(Tmp,2);
                               _Mg_Mod_Pow_LInt(a,tmp,modulo,MgStruct,Output);
                               end
else begin
                        {The Modulo is in the form 4k+1 so we use the Tonelli-Shanks Algorithm}
     _Sub_LInt(modulo,1,MinusOne);
     _HCopy_LInt(MinusOne,MinusOneDiv2);
     _Shr_LInt(Minusonediv2,1);
     _HCopy_LInt(2,g);
     _Mg_Mod_Pow_LInt(g,Minusonediv2,Modulo,MgStruct,tmp);
     while _Compare_LInt(tmp,MinusOne)<>0 do begin
                                           _Inc_LInt(g,1);
                                           _Mg_Mod_Pow_LInt(g,Minusonediv2,Modulo,MgStruct,tmp);
                                           end;
     _HCopy_LInt(MinusOneDiv2,e1);
     e2.data.i32[-1]:=0;
     While(e1.data.i32[0] and 1=0) do begin
                                 _Shr_LInt(e1,1);
                                 _Shr_LInt(e2,1);
                                 _Mg_Mod_Pow_LInt(a,e1,Modulo,MgStruct,tmp);
                                 _Mg_Mod_Pow_LInt(g,e2,Modulo,MgStruct,tmp2);
                                 _Mul_LInt(tmp,tmp2,tmp3);
                                 _Mod_LInt(tmp3,Modulo,tmp);
                                 if _Compare_LInt(tmp,MinusOne)=0 then _Add_LInt(e2,MinusOneDiv2,e2)
                                 end;
     _Inc_LInt(e1,1);
     _Shr_LInt(e1,1);
     _Mg_Mod_Pow_LInt(a,e1,Modulo,MgStruct,tmp2);
     _Shr_LInt(e2,1);
     _Mg_Mod_Pow_LInt(g,e2,Modulo,MgStruct,tmp3);
     _Mul_LInt(tmp2,tmp3,tmp);
     _Mod_LInt(tmp,Modulo,Output);
     end;
end;

{***************************************************************************************}
{              Convert a Big integer to a decimal string                                }
{***************************************************************************************}
Function _LInt_To_Str(a:LInt):String;
var     s,sig:string;
        rem,bilion:LInt;
        SaveSig:Word;
begin
bilion.data.i32[-1]:=1;
bilion.data.i32[0]:=1000000000;
result:='';
SaveSig:=a.data.i16[-1];
if SaveSig<>0 then sig:='-' else sig:='';
a.data.i16[-1]:=0;
repeat
     _Div_Mod_LInt(a,bilion,a,rem);
     if rem.data.i16[-2]<>0 then s:=inttostr(rem.data.i32[0])
     else s:='0';
     while length(s)<9 do s:='0'+s;
     Result:=s+Result;
   until(a.data.i32[-1]and $ffff=1)or(a.data.i32[-1]and $ffff=0)and(a.data.i32[0]<=999999999);
Result:=inttostr(a.data.i32[0])+result;
while (result[1]='0')and(length(Result)>1) do result:=copy(result,2,length(result));
result:=sig+result;
a.data.i16[-1]:=SaveSig;
end;

{***************************************************************************************}
{              Convert a Big integer to a Hexadecimal string                            }
{***************************************************************************************}
Function _LInt_To_Hex(a:LInt):String;
var i:integer;
    s:string;
begin
Result:='';
if a.data.i16[-2] =0 then result:='0'
else
for i:=0 to (a.data.i16[-2]-1) mod (MaxLen shl 5) do begin
                           s:=inttohex(a.data.i32[i],8);
                           Result:=s+Result;
                           end;
while (Result[1]='0')and(length(Result)>1) do Result:=copy(Result,2,length(Result));
Result:='0x'+Result;
if _IsNeg(a) then result:='-'+Result;
end;

{***************************************************************************************}
{              Convert a ByteArray to a Hexadecimal string                              }
{***************************************************************************************}
Function ArrayToHex(a:Tbytes):wideString;
var i:integer;
    s:widestring;
begin
Result:='';
if Length(a)=0 then result:='0'
else
for i:=0 to length(a)-1 do begin
                           s:=inttohex(a[i],2);
                           Result:=s+Result;
                           end;
result:='0x'+result;
end;
{***************************************************************************************}
{              Convert a Hexadecimal string to ByteArray                                }
{***************************************************************************************}
Function HexToArray(a:widestring):Tbytes;
var i:integer;
    s:string;
begin
s:=a;
if (s[1]='$') then s:=copy(S,2,length(s))
else if (s[1]='0')and(s[2]='x') then s:=copy(S,3,length(s));
setlength(result,(length(s) div 2+length(s) mod 2));
i:=0;
While (s<>'') do begin
                 if length(s)>2 then begin
                                     Result[i]:=strtoint('$'+copy(s,length(s)-1,2));
                                     s:=copy(s,1,length(s)-2);
                                     inc(i);
                                     end
                 else begin
                      Result[i]:=strtoint('$'+s);
                      s:='';
                      end;
                 end;
end;


{***************************************************************************************}
{              Convert Bytes Array integer to a Big Integer                             }
{***************************************************************************************}
Function ByteArrayToLint(a:TBytes):LInt;
var i,Size:integer;
begin
if Length(a) mod 4 =0 then Size:=Length(a) div 4 else Size:=Length(a) div 4 +1;
Result.Data.i16[-2]:=Size;
Result.Data.i16[-1]:=0;
Result.Data.i32[Size-1]:=0;
for i:=0 to Length(a)-1 do Result.Data.i8[i]:=a[i];
end;

{***************************************************************************************}
{              Convert a Big Integer  to a Bytes Array integer                          }
{***************************************************************************************}
Function LIntToByteArray(a:LInt):Tbytes;
var i,Size:integer;
begin
if a.BitLength mod 8 =0 then Size:=a.BitLength div 8 else Size:=a.BitLength div 8+1;
SetLength(Result,Size);
for i:=0 to Size-1 do Result[i]:=a.Data.i8[i];
end;

{***************************************************************************************}
{              Convert a Hexadecimal string to a Bige integer                           }
{***************************************************************************************}
Function _Hex_To_LInt(a:String):LInt;
var i:integer;
    s:string;
    signe:Dword;
begin
signe:=0;
s:=a;
if s[1]='-' then begin
                 s:=copy(S,2,length(s));
                 signe:=$ffff0000;
                 end;
if (s[1]='$') then s:=copy(S,2,length(s))
else if (s[1]='0')and(s[2]='x') then s:=copy(S,3,length(s));
i:=0;
While (s<>'') do begin
                 if length(s)>8 then begin
                                     Result.data.i32[i]:=strtoint('$'+copy(s,length(s)-7,8));
                                     s:=copy(s,1,length(s)-8);
                                     inc(i);
                                     end
                 else begin
                      Result.data.i32[i]:=strtoint('$'+s);
                      s:='';
                      end;
                 end;
Result.data.i32[-1]:=signe or (i+1);
if (Result.data.i16[-2] =1)and(Result.Data.i32[0]=0) then Result.Data.i16[-2]:=0;

end;
{***************************************************************************************}
{              Convert a decimal string to a Big integer                                }
{***************************************************************************************}
Function _Str_To_LInt(s:string):LInt;
var str:string;
    signe:Dword;
    million,tmp,t,z,h:LInt;
begin
signe:=0;
if s[1]='-' then begin
                 s:=copy(S,2,length(s));
                 signe:=$ffff0000;
                 end;
Million.data.i32[-1]:=1;
t.data.i32[-1]:=1;
tmp.data.i32[-1]:=1;
t.data.i32[0]:=1;
Million.data.i32[0]:=1000000;
if length(s)<=7 then begin
                     Result.data.i32[0]:=strtoint(s);
                     if Result.data.i32[0]<>0 then Result.data.i32[-1]:=1
                     else Result.data.i32[-1]:=0;
                     end
else begin
     str:=s;
     Result.data.i32[-1]:=1;
     Result.data.i32[0]:=0;
        repeat
           tmp.data.i32[0]:=strtoint(copy(str,length(str)-5,length(str)));
           str:=copy(str,1,length(str)-6);
           _Mul_LInt(tmp,t,z);
           _Add_LInt(Result,z,Result);
           _Mul_LInt(t,million,h);
           _HCopy_LInt(h,t);
         until(str='');
     end;
Result.data.i32[-1]:=Result.data.i32[-1] or signe;
end;

{***************************************************************************************}
Function IsLIntOdd(a:LInt):boolean;
begin
if not _isnull(a) then result:=(a.Data.i32[0] and 1)>0
else Result:=true;
end;

{***************************************************************************************}
Function IsLIntEven(a:LInt):boolean;
begin
result:=(a.Data.i32[0] and 1)=0;
end;

{***************************************************************************************}
function LIntToNAF(val:LInt):LIntArrayForm;
var i:integer;
    a:Lint;
begin
a:=val.Absolute;
Setlength(result,a.BitLength+1);
i:=0;
while (a<>0) do begin
               if a.isOdd then begin
                               Result[i]:=2-(a.Data.i32[0] and 3);
                               a:=a-Result[i];
                               end
               else Result[i]:=0;
               _Shr_LInt(a,1);
               inc(i);
               end;
if Result[length(Result)-1]=0 then setlength(Result,Length(Result)-1);
end;

{***************************************************************************************}
function LIntToIntArray(a:LInt):LIntArrayForm;
var i:integer;
begin
Setlength(result,a.BitLength);
i:=0;
while (a<>0) do begin
               if a.isOdd then Result[i]:=1
               else Result[i]:=0;
               _Shr_LInt(a,1);
               inc(i);
               end;
end;

{***************************************************************************************}
function ArrayHammingWeight(L:LIntArrayForm):Dword;
var i:integer;
begin
Result:=0;
for i:=0 to length(L)-1 do if L[i]<>0 then inc(Result);
end;

{***************************************************************************************}
function HammingWeight(L:Lint;Negative:Boolean=False):Dword;
var t1:LIntArrayForm;
    i:integer;
begin
if Negative then t1:=LIntToNAF(L.Absolute)
else t1:=LIntToIntArray(L.Absolute);
Result:=0;
for i:=0 to length(t1)-1 do if t1[i]<>0 then inc(Result);
end;

{***************************************************************************************}
function NAFToLInt(naf:LIntArrayForm):LInt;
var i:integer;
begin
Result:=0;
for i:=length(naf)-1 downto 0 do begin
                             if naf[i]=1 then Result:=Result+(LInt(1) shl i)
                             else if naf[i]=-1 then Result:=Result-(LInt(1) shl i);
                             end;
end;

{***************************************************************************************}
{              Compute the GCD of two Big integers using Binary Algo                    }
{***************************************************************************************}
procedure BinaryGCDLInt(a,b:LInt;var Result:LInt);inline;
var u,v,g,t:LInt;
begin
Result.data.i32[-1]:=0;
if (_IsNull(a))or (_IsNull(b)) then exit
else begin
     _HCopy_LInt(a,u);
     _HCopy_LInt(b,v);
     _HCopy_LInt(1,g);
     While(u.data.i32[0] and 1=0)and(v.data.i32[0] and 1=0) do begin
                                                     _Shr_LInt(u,1);
                                                     _Shr_LInt(v,1);
                                                     _Shl_LInt(g,1);
                                                     end;
     While not(_isNull(v)) do begin
                             if (u.data.i32[0] and 1=0) then _Shr_LInt(u,1)
                             else if (v.data.i32[0] and 1=0) then _Shr_LInt(v,1)
                             else begin
                                  _Sub_LInt(u,v,t);
                                  t.data.i16[-1]:=0;
                                  _Shr_LInt(t,1);
                                  if _Compare_LInt(u,v)<=0 then _HCopy_LInt(t,v) else _HCopy_LInt(t,u);
                                  end;
                             end;
     _Mul_LInt(u,g,Result);
     end;
end;

{**********************************************************************************************************************}
function LInt.ToDecimalString:string;
begin
Result:=_LInt_To_Str(Self);
end;

{**********************************************************************************************************************}
function LInt.ToHexString:String;
begin
Result:=_LInt_To_Hex(Self);
end;

{**********************************************************************************************************************}
function LInt.ToNafArray: LIntArrayForm;
begin
Result:=LIntToNAF(Self);
end;

{**********************************************************************************************************************}
function LInt.ToIntArray: LIntArrayForm;
begin
Result:=LIntToIntArray(Self);
end;


{**********************************************************************************************************************}
Function LInt.InversModulo(Modulo:LInt):LInt;
begin
_Inv_Mod_LInt(Self,Modulo,Result);
end;

{**********************************************************************************************************************}
function LInt.IsEven: boolean;
begin
result:=IsLIntEven(self);
end;

{**********************************************************************************************************************}
function LInt.IsZero: boolean;
begin
result:=_IsNull(self);
end;

{**********************************************************************************************************************}
function LInt.IsOdd: boolean;
begin
result:=IsLIntOdd(Self);
end;

{**********************************************************************************************************************}
function LInt.IsASqrMod(p: LInt;MgStruct:Pointer): boolean; /// p is supposed to be a prime number
var g: LInt;
begin
_Mg_Mod_Pow_LInt(self,((p-1) shr 1),p,MgStruct,g);
if (_Compare_LInt(Self mod p,0)<>0)and(_Compare_LInt(g,1)<>0) then Result:=false
else result:=true;
end;

{**********************************************************************************************************************}
procedure LInt.SetToRandom(Limit: LInt;Seed:Word=0);
begin
GetRandomLIntLowerThan(Self,Limit,seed);
end;

{**********************************************************************************************************************}
procedure LInt.SetToRandomOnBits(nbits: Integer);
begin
GetRandomLIntOnBits(Self,nbits);
end;

{**********************************************************************************************************************}
Function LInt.Sqr:LInt;
begin
_Sqr_LInt(Self,Result);
end;

{**********************************************************************************************************************}
function LInt.Sqrt(p:LInt;MgStruct:Pointer): LInt;
begin
if not ModSquareLInt(Self,p, result,PtrMontgomeryStruct(Mgstruct)) then raise ERangeError.Create( 'This number is not a Quadratic Residue......');
end;

{**********************************************************************************************************************}
class operator LInt.Subtract(const Left: LInt; Right: Integer): LInt;
var x:LInt;
begin
x.Data.i16[-2]:=1;
if right>=0 then x.Data.i16[-1]:=0 else x.Data.i16[-1]:=$ffff;
x.Data.i32[0]:=abs(Right);
_Sub_LInt(Left,x,Result);
end;

{**********************************************************************************************************************}
function LInt.BitLength:word;
begin
result:=_Num_Of_Bits_LInt(Self);
end;

{**********************************************************************************************************************}
class operator LInt.Divide(const Left: LInt; Right: LInt): LInt;
var tmp:LInt;
begin
_Div_Mod_LInt(Left,Right,Result,tmp);
end;

{**********************************************************************************************************************}
function LInt.TestBit(i:integer):boolean;
begin
result:=_Is_BitSet_At(Self,i);
end;

{**********************************************************************************************************************}
class operator LInt.Add(const Left, Right: LInt): LInt;
begin
_Add_LInt(Left,Right,Result);
end;

{**********************************************************************************************************************}
class operator LInt.Add(const Left:LInt; Right: word): LInt;
begin
_HCopy_LInt(Left,Result);
_Inc_LInt(Result,Right);
end;

{**********************************************************************************************************************}
function LInt.Absolute: LInt;
begin
Result:=Self;
Result.Data.i16[-1]:=0;
end;

{**********************************************************************************************************************}
class operator LInt.Add(const Left:Word; Right: LInt): LInt;
begin
_HCopy_LInt(Right,Result);
_Inc_LInt(Result,Left);
end;

{**********************************************************************************************************************}
class operator LInt.Subtract(const Left, Right: LInt): LInt;
begin
_Sub_LInt(Left,Right,Result);
end;

{**********************************************************************************************************************}
class operator LInt.Multiply(const Left, Right: LInt): LInt;
begin
_Mul_LInt(Left,Right,Result);
end;

{**********************************************************************************************************************}
class operator LInt.Multiply(const Left: LInt; Right: Word): LInt;
var tmp:LInt;
begin
tmp:=0;
_Inc_LInt(tmp,Right);
_Mul_LInt(tmp,Left,Result);
end;

{**********************************************************************************************************************}
class operator LInt.Multiply(Left: integer; const Right: LInt): LInt;
var tmp:LInt;
begin
tmp.Data.i16[-2]:=1;
if left>0 then tmp.Data.i16[-1]:=0 else tmp.Data.i16[-1]:=$ffff;
tmp.Data.i32[0]:=abs(Left);
_Mul_LInt(tmp,Right,Result);
end;

{**********************************************************************************************************************}
class operator LInt.Implicit(const Value: string): LInt;
var neg:boolean;
    S:string;
begin
S:=Value;
neg:=S[1]='-';
if neg then s:=copy(s,2,length(s));
if (S[1]='$') then Result:=_Hex_To_LInt(Copy(S,2,length(S)))
else if(S[1]='0')and(S[2]='x') then Result:=_Hex_To_LInt(Copy(S,3,length(S)))
else Result:=_Str_To_LInt(S);
if neg then Result.Data.i16[-1]:=$ffff;
end;

{**********************************************************************************************************************}
class operator LInt.Implicit(const Value: Longint): LInt;
begin
if Value=0 then Result.Data.i32[-1]:=0
else begin
     Result.Data.i16[-2]:=1;
     if value>=0 then  Result.Data.i16[-1]:=0 else Result.Data.i16[-1]:=$ffff;
     Result.Data.i32[0]:=abs(Value);
     end;
end;


{**********************************************************************************************************************}
class operator LInt.Modulus(const Left: LInt; Right: LInt): LInt;
begin
_Mod_LInt(Left,Right,Result);
end;

{**********************************************************************************************************************}
class operator LInt.LeftShift(const Value: LInt; Shift: Integer): LInt;
begin
_HCopy_LInt(Value,Result);
_Shl_LInt(Result,Shift);
end;

{**********************************************************************************************************************}
class operator LInt.RightShift(const Value: LInt; Shift: Integer): LInt;
begin
_HCopy_LInt(Value,Result);
_Shr_LInt(Result,Shift);
end;

{**********************************************************************************************************************}
class operator LInt.Equal(const Left, Right: LInt): Boolean;
begin
Result:=(_Compare_LInt(Left,Right)=0);
end;

{**********************************************************************************************************************}
class operator LInt.Equal(const Left:LInt; Right: DWord): Boolean;
var tmp:Dword;
begin
if Right=0 then Result:=_IsNull(Left)
else begin
     tmp:=Left.Data.i32[0];
     if Left.Data.i16[-1]<>0 then tmp:=-tmp;
     Result:=(Left.Data.i16[-2]=1)and(tmp=Right);
     end;
end;

{**********************************************************************************************************************}
class operator LInt.NotEqual(const Left, Right: LInt): Boolean;
begin
Result:=(_Compare_LInt(Left,Right)<>0);
end;

{**********************************************************************************************************************}
class operator LInt.NotEqual(const Left:LInt; Right: Dword): Boolean;
begin
if Right=0 then Result:=not _IsNull(Left)
else Result:=(Left.Data.i16[-2]=1)and(Left.Data.i32[0]<>Right);
end;

{**********************************************************************************************************************}
function LInt.Pow(e: longword): LInt;
begin
_Pow_LInt(Self,Result,e);
end;

{**********************************************************************************************************************}
function LInt.PowMod(e, p: LInt;MgStr:Pointer=nil ): LInt;
begin
if MgStr=nil then _Mod_Pow_LInt(Self,e,p,Result)
else _Mg_Mod_Pow_LInt(Self,e,p, PtrMontgomeryStruct(Mgstr),Result);
end;

{**********************************************************************************************************************}
class operator LInt.GreaterThan(const Left, Right: LInt): Boolean;
begin
Result:=(_Compare_LInt(Left,Right)=1);
end;

{**********************************************************************************************************************}
class operator LInt.GreaterThan(const Left:LInt; Right: Dword): Boolean;
begin
if _IsNeg(Left) then Result:=false
else if Right=0 then Result:=not _IsNull(Left)
else Result:=((Left.Data.i16[-2]=1)and(Left.Data.i32[0]>Right)) or(Left.Data.i16[-2]>1);
end;

{**********************************************************************************************************************}
class operator LInt.GreaterThanOrEqual(const Left, Right: LInt): Boolean;
begin
Result:=(_Compare_LInt(Left,Right)>=0);
end;

{**********************************************************************************************************************}
class operator LInt.GreaterThanOrEqual(const Left:LInt; Right: Dword): Boolean;
begin
if _IsNeg(Left) then result:=false
else Result:=((Left.Data.i16[-2]=1)and(Left.Data.i32[0]>=Right))or(Left.Data.i16[-2]>1);
end;

{**********************************************************************************************************************}
class operator LInt.LessThan(const Left, Right: LInt): Boolean;
begin
Result:=(_Compare_LInt(Left,Right)<0);
end;

{**********************************************************************************************************************}
class operator LInt.LessThan(const Left:LInt; Right: Dword): Boolean;
begin
if Right=0 then Result:=_IsNeg(Left)
else Result:=(Left.Data.i16[-2]=1)and(Left.Data.i32[0]<Right);
end;

{**********************************************************************************************************************}
class operator LInt.LessThanOrEqual(const Left:LInt; Right: Dword): Boolean;
begin
if Right=0 then Result:=_IsNeg(Left)or _IsNull(Left)
else Result:=(Left.Data.i16[-2]=1)and(Left.Data.i32[0]<=Right);
end;

{**********************************************************************************************************************}
class operator LInt.LessThanOrEqual(const Left, Right: LInt): Boolean;
begin
Result:=(_Compare_LInt(Left,Right)<=0);
end;

{**********************************************************************************************************************}
class operator LInt.Negative(const Value: LInt): LInt;
begin
Result:=value;
if Result.Data.i16[-1]=0 then Result.Data.i16[-1]:=$ffff else Result.Data.i16[-1]:=0;
end;

end.



