#include "DxLib.h"
#include "const.h"
#include "calc.h"
#include <math.h>
#include <stdlib.h>

#define PAGE_MAX 6

struct planet data[P_MAX];


double wt = 0.0;
int viewtype = -1;
float objsize = 4.0f;
double dt = UNIT_MIN*UNIT_HOUR/4;
int metn = 0;
int draw=16;
double eps=0.000001;

int p_main() {
	double r[P_MAX][6];
	int i,j;
	
	for(i=0;i<pn;i++) {
		if(metn==1) {
			Euler(i,dt,r[i]);
		} else if(metn== 2) {
			impEuler(i,dt,r[i]);
		} else if(metn== 3) {
			MidPoint(i,dt,r[i]);
		} else if(metn== 4) {
			Ralston(i,dt,r[i]);
		} else if(metn== 5) {
			Kutta3o(i,dt,r[i]);
		} else if(metn== 6) {
			RungeKutta4_38(i,dt,r[i]);
		} else if(metn== 7) {
			RungeKuttaFehlberg(i,dt,eps,r[i]);
		} else if(metn== 8) {
			VelocityVerlet(i,dt,r[i]);
		} else if(metn== 9) {
			EulerRichardson(i,dt,r[i]);
		} else if(metn==10) {
			Leapfrog(i,dt,r[i]);
		} else {
			RungeKutta(i,dt,r[i]);
		}
	}
	for(i=0;i<pn;i++) {
		for(j=0;j<3;j++) {
			data[i].x[j] = r[i][j  ];
			data[i].v[j] = r[i][j+3];
		}
	}
	
	if(viewtype >= 0 && viewtype < pn)
		center(viewtype);
		
	
	wt += dt;
	
	return 0;
}

// (x,y)の点を(mx,my)を中心にang角回転する
void rotate(float *x, float *y, const float ang, const float mx, const float my){
    const float ox = *x-mx, oy = *y-my;
    *x =  ox * cos(ang) + oy * sin(ang);
    *y = -ox * sin(ang) + oy * cos(ang);
    *x += mx;
    *y += my;
}

void print_planets(float *x,double wt,float *rot,float *move,VECTOR CameraPos,int rn,int mn,char hide) {
	int i,j;
	// 画面をクリア
	ClearDrawScreen() ;
	
	// ３Ｄ空間上に球を描画する
	for(i=0;i<pn;i++) {
		for(j=0;j<3;j++)
			x[j]=(float)length(data[i].x[j]);
		DrawSphere3D( VGet( x[0], x[1], x[2] ), (float)radius(x[0]-CameraPos.x,x[1]-CameraPos.y,CameraPos.z-x[2]) / 200, 32, GetColor( data[i].c[0],data[i].c[1],data[i].c[2] ), GetColor( 255, 255, 255 ), TRUE ) ;
		DrawSphere3D( VGet( x[0], x[1], x[2] ), objsize, 32, GetColor( data[i].c[0],data[i].c[1],data[i].c[2] ), GetColor( 255, 255, 255 ), TRUE ) ;
	}
	if(hide<2)
		DrawFormatString( 0, 0, GetColor( 0, 255, 0 ), "time   : %8.3f [day]", time(wt) ); // 文字を描画する
	if(hide<1) {
		DrawFormatString( 0,20, GetColor( 0, 255, 0 ), "rotate : %8.3f [rad]", rot[rn] );
		DrawFormatString( 0,40, GetColor( 0, 255, 0 ), "move   : %8.3f [dot]", move[mn] );
	}
}

// プログラムは WinMain から始まります
int WINAPI WinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance,
						LPSTR lpCmdLine, int nCmdShow )
{
	char key[256];
	char buf[32];
	VECTOR CameraPos ;
	VECTOR ViewPoint ;
	VECTOR UPVector  ;
	VECTOR DirectX   ;
	VECTOR DirectY   ;
	VECTOR DirectZ   ;
	int mn=2,rn=2,i,j;
	int count,wi,he;
	unsigned int c=0;
	char hide = 1;
	double r;
	float x[3];
	int page=0;
	int tmpkey[256];
	float move[]={1000.0f,100.0f,10.0f,1.0f,0.1f,0.01f};
	float rot[] ={(float)DX_PI_F/15,(float)DX_PI_F/36,(float)DX_PI_F/90,(float)DX_PI_F/180};
	char method[][64]={"4-stage Runge-Kutta method","Euler method","improved Euler method","midpoint method","Ralston's method","Kutta's 3rd-order method","3/8-rule 4th-order method","Runge-Kutta-Fehlberg method","Velocity Verlet method","Euler Richardson method","Leap-frog scheme"};
	
	// タイトルを test に変更
	SetMainWindowText( "Planet (DirectX version)" ) ;
	
	// ウインドウモードに変更
	ChangeWindowMode(TRUE);
	if( DxLib_Init() == -1 )		// ＤＸライブラリ初期化処理
	{
		return -1 ;			// エラーが起きたら直ちに終了
	}
	// 画面モードの設定
	//SetGraphMode(960,720,32);
	
	//GetWindowSize(&wi,&he);
	//SetDrawArea(0,0,he+1,he+1);
	//SetCameraScreenCenter(he/2,he/2);
	
	input();
	center(viewtype);
        
	// 描画先を裏画面にする
	SetDrawScreen( DX_SCREEN_BACK ) ;
	
	// ライティングの計算をしないように設定を変更
	SetUseLighting( FALSE ) ;

        // Ｚバッファを有効にする
	SetUseZBuffer3D( TRUE ) ;

	// Ｚバッファへの書き込みを有効にする
	SetWriteZBuffer3D( TRUE ) ;
    
	// カメラの座標を初期化
	CameraPos.x = 0.0f ;
	CameraPos.y = 0.0f ;
	CameraPos.z = -500.0f;
	ViewPoint.x = 0.0f;
	ViewPoint.y = 0.0f;
	ViewPoint.z = 0.0f;
	UPVector.x  = 0.0f;
	UPVector.y  = 1.0f;
	UPVector.z  = 0.0f;
	
	
        // ＥＳＣキーが押されるかウインドウが閉じられるまでループ
	while( ProcessMessage() == 0 && CheckHitKey( KEY_INPUT_ESCAPE ) == 0 )
	{
		p_main();
		
		if(c!=0){
			c=(c+1)%draw;
			continue;
		}
		
		print_planets(x,wt,rot,move,CameraPos,rn,mn,hide);
		
		// 画面をクリア
		//ClearDrawScreen() ;
		
		// ３Ｄ空間上に球を描画する
		//for(i=0;i<pn;i++) {
		//	for(j=0;j<3;j++)
		//		x[j]=(float)length(data[i].x[j]);
		//	DrawSphere3D( VGet( x[0], x[1], x[2] ), (float)radius(x[0]-CameraPos.x,x[1]-CameraPos.y,CameraPos.z-x[2]) / 200, 32, GetColor( data[i].c[0],data[i].c[1],data[i].c[2] ), GetColor( 255, 255, 255 ), TRUE ) ;
		//	DrawSphere3D( VGet( x[0], x[1], x[2] ), objsize, 32, GetColor( data[i].c[0],data[i].c[1],data[i].c[2] ), GetColor( 255, 255, 255 ), TRUE ) ;
		//}
		//if(hide<2)
		//	DrawFormatString( 0, 0, GetColor( 0, 255, 0 ), "time   : %8.3f [day]", time(wt) ); // 文字を描画する
		//if(hide<1) {
		//	DrawFormatString( 0,20, GetColor( 0, 255, 0 ), "rotate : %8.3f [rad]", rot[rn] );
		//	DrawFormatString( 0,40, GetColor( 0, 255, 0 ), "move   : %8.3f [dot]", move[mn] );
		//}
		
		//DrawFormatString( 0,60, GetColor( 0, 255, 0 ), "(%f,%f,%f)", ViewPoint.x, ViewPoint.y, ViewPoint.z ); // 文字を描画する
		
		GetHitKeyStateAll( key );
		
		if(( key[KEY_INPUT_LSHIFT] | key[KEY_INPUT_RSHIFT] ) == 1) {
			if((key[KEY_INPUT_LCONTROL] | key[KEY_INPUT_LCONTROL]) == 1) {
				if( key[KEY_INPUT_LEFT ] == 1 )
				{
					rotate(&CameraPos.x,&CameraPos.z, rot[rn],ViewPoint.x,ViewPoint.z);
				}
				if( key[KEY_INPUT_RIGHT ] == 1 )
				{
					rotate(&CameraPos.x,&CameraPos.z,-rot[rn],ViewPoint.x,ViewPoint.z);
				}
				if( key[KEY_INPUT_UP ] == 1 )
				{
					rotate(&CameraPos.y,&CameraPos.z,-rot[rn],ViewPoint.y,ViewPoint.z);
					rotate(&UPVector.y ,&UPVector.z ,-rot[rn],0.0f,0.0f);
				}
				if( key[KEY_INPUT_DOWN ] == 1 )
				{
					rotate(&CameraPos.y,&CameraPos.z, rot[rn],ViewPoint.y,ViewPoint.z);
					rotate(&UPVector.y ,&UPVector.z , rot[rn],0.0f,0.0f);
				}
				if( key[KEY_INPUT_A ] == 1 )
				{
					rotate(&CameraPos.x,&CameraPos.y,-rot[rn],ViewPoint.x,ViewPoint.y);
					rotate(&UPVector.x ,&UPVector.y ,-rot[rn],0.0f,0.0f);
				}
				if( key[KEY_INPUT_D ] == 1 )
				{
					rotate(&CameraPos.x,&CameraPos.y, rot[rn],ViewPoint.x,ViewPoint.y);
					rotate(&UPVector.x ,&UPVector.y , rot[rn],0.0f,0.0f);
				}
				for(i=0;i<256;i++) {
					if(key[i] == 1) {
						tmpkey[i]++;
					} else {
						tmpkey[i]=0;
					}
				}
				if( tmpkey[KEY_INPUT_Z ] == 1 )
				{
					rn=(rn+1)%(sizeof(rot)/sizeof(rot[0]));
				}
				if( tmpkey[KEY_INPUT_X ] == 1 )
				{
					rn=(rn-1+(sizeof(rot)/sizeof(rot[0])))%(sizeof(rot)/sizeof(rot[0]));
				}
			} else {
				if( key[KEY_INPUT_LEFT ] == 1 )
				{
					rotate(&ViewPoint.x,&ViewPoint.z,-rot[rn],CameraPos.x,CameraPos.z);
				}
				if( key[KEY_INPUT_RIGHT ] == 1 )
				{
					rotate(&ViewPoint.x,&ViewPoint.z, rot[rn],CameraPos.x,CameraPos.z);
				}
				if( key[KEY_INPUT_UP ] == 1 )
				{
					rotate(&ViewPoint.y,&ViewPoint.z, rot[rn],CameraPos.y,CameraPos.z);
					rotate(&UPVector.y ,&UPVector.z , rot[rn],0.0f,0.0f);
				}
				if( key[KEY_INPUT_DOWN ] == 1 )
				{
					rotate(&ViewPoint.y,&ViewPoint.z,-rot[rn],CameraPos.y,CameraPos.z);
					rotate(&UPVector.y ,&UPVector.z ,-rot[rn],0.0f,0.0f);
				}
				if( key[KEY_INPUT_A ] == 1 )
				{
					rotate(&ViewPoint.x,&ViewPoint.y,-rot[rn],CameraPos.x,CameraPos.y);
					rotate(&UPVector.x ,&UPVector.y ,-rot[rn],0.0f,0.0f);
				}
				if( key[KEY_INPUT_D ] == 1 )
				{
					rotate(&ViewPoint.x,&ViewPoint.y, rot[rn],CameraPos.x,CameraPos.y);
					rotate(&UPVector.x ,&UPVector.y , rot[rn],0.0f,0.0f);
				}
				for(i=0;i<256;i++) {
					if(key[i] == 1) {
						tmpkey[i]++;
					} else {
						tmpkey[i]=0;
					}
				}
				if( tmpkey[KEY_INPUT_Z ] == 1 || (tmpkey[KEY_INPUT_Z ] >= 60 && tmpkey[KEY_INPUT_Z ] % 60 == 0) )
				{
					rn=(rn+1)%(sizeof(rot)/sizeof(rot[0]));
				}
				if( tmpkey[KEY_INPUT_X ] == 1 || (tmpkey[KEY_INPUT_X ] >= 60 && tmpkey[KEY_INPUT_X ] % 60 == 0) )
				{
					rn=(rn-1+(sizeof(rot)/sizeof(rot[0])))%(sizeof(rot)/sizeof(rot[0]));
				}
				if( tmpkey[KEY_INPUT_H ] == 1 )
				{
					hide=(hide+1)%2;
				}
			}
		} else if((key[KEY_INPUT_LALT] | key[KEY_INPUT_RALT]) == 1) {
			DirectY.x = UPVector.x;
			DirectY.y = UPVector.y;
			DirectY.z = UPVector.z;
			r = radius(ViewPoint.x-CameraPos.x,ViewPoint.y-CameraPos.y,CameraPos.z-ViewPoint.z);
			DirectZ.x = (ViewPoint.x-CameraPos.x)/r;
			DirectZ.y = (ViewPoint.y-CameraPos.y)/r;
			DirectZ.z = (ViewPoint.z-CameraPos.z)/r;
			DirectX.x = DirectY.y*DirectZ.z-DirectZ.y*DirectY.z;
			DirectX.y = DirectY.z*DirectZ.x-DirectZ.z*DirectY.x;
			DirectX.z = DirectY.x*DirectZ.y-DirectZ.y*DirectY.x;
			if( key[KEY_INPUT_LEFT ] == 1 )
			{
				ViewPoint.x += (float)move[mn] * DirectX.x ;
				ViewPoint.y += (float)move[mn] * DirectX.y ;
				ViewPoint.z += (float)move[mn] * DirectX.z ;
			}
			if( key[KEY_INPUT_RIGHT ] == 1 )
			{
				ViewPoint.x -= (float)move[mn] * DirectX.x ;
				ViewPoint.y -= (float)move[mn] * DirectX.y ;
				ViewPoint.z -= (float)move[mn] * DirectX.z ;
			}
			if( key[KEY_INPUT_UP ] == 1 )
			{
				ViewPoint.x += (float)move[mn] * DirectY.x ;
				ViewPoint.y += (float)move[mn] * DirectY.y ;
				ViewPoint.z += (float)move[mn] * DirectY.z ;
			}
			if( key[KEY_INPUT_DOWN ] == 1 )
			{
				ViewPoint.x -= (float)move[mn] * DirectY.x ;
				ViewPoint.y -= (float)move[mn] * DirectY.y ;
				ViewPoint.z -= (float)move[mn] * DirectY.z ;
			}
			if( key[KEY_INPUT_W ] == 1 )
			{
				ViewPoint.x += (float)move[mn] * DirectZ.x ;
				ViewPoint.y += (float)move[mn] * DirectZ.y ;
				ViewPoint.z += (float)move[mn] * DirectZ.z ;
			}
			if( key[KEY_INPUT_S ] == 1 )
			{
				ViewPoint.x -= (float)move[mn] * DirectZ.x ;
				ViewPoint.y -= (float)move[mn] * DirectZ.y ;
				ViewPoint.z -= (float)move[mn] * DirectZ.z ;
			}
			for(i=0;i<256;i++) {
				if(key[i] == 1) {
					tmpkey[i]++;
				} else {
					tmpkey[i]=0;
				}
			}
			if( tmpkey[KEY_INPUT_Z ] == 1 || (tmpkey[KEY_INPUT_Z ] >= 60 && tmpkey[KEY_INPUT_Z ] % 60 == 0) )
			{
				mn=(mn+1)%(sizeof(move)/sizeof(move[0]));
			}
			if( tmpkey[KEY_INPUT_X ] == 1 || (tmpkey[KEY_INPUT_X ] >= 60 && tmpkey[KEY_INPUT_X ] % 60 == 0) )
			{
				mn=(mn-1+(sizeof(move)/sizeof(move[0])))%(sizeof(move)/sizeof(move[0]));
			}
		} else {
			DirectY.x = UPVector.x;
			DirectY.y = UPVector.y;
			DirectY.z = UPVector.z;
			r = radius(ViewPoint.x-CameraPos.x,ViewPoint.y-CameraPos.y,CameraPos.z-ViewPoint.z);
			DirectZ.x = (ViewPoint.x-CameraPos.x)/r;
			DirectZ.y = (ViewPoint.y-CameraPos.y)/r;
			DirectZ.z = (ViewPoint.z-CameraPos.z)/r;
			DirectX.x = DirectY.y*DirectZ.z-DirectZ.y*DirectY.z;
			DirectX.y = DirectY.z*DirectZ.x-DirectZ.z*DirectY.x;
			DirectX.z = DirectY.x*DirectZ.y-DirectZ.y*DirectY.x;
			
			if( key[KEY_INPUT_C ] == 1 ) {
				buf[0]=0;
				while( ProcessMessage() == 0 && CheckHitKey( KEY_INPUT_ESCAPE ) == 0 ) {
					// 画面をクリア
					ClearDrawScreen() ;
					
					if(page == 0) {
						DrawFormatString( 0, 0, GetColor( 0, 255, 0 ), "The center object of this window : %s", buf );
						DrawFormatString( 0,20, GetColor( 0, 255, 0 ), "(The center of gravity : -1)" );
						DrawFormatString(20,40, GetColor( 0, 255, 0 ), "Now : %d", viewtype );
					} else if(page == 1) {
						DrawFormatString( 0, 0, GetColor( 0, 255, 0 ), "dt[s] : %s", buf );
						DrawFormatString(20,20, GetColor( 0, 255, 0 ), "Now : %f", dt );
					} else if(page == 2) {
						DrawFormatString( 0, 0, GetColor( 0, 255, 0 ), "Method : %s", buf );
						for(i=0;i<sizeof(method)/sizeof(method[0]);i++)
							DrawFormatString( 0,20*(i+1), GetColor( 0, 255, 0 ), "%d : %s", i, method[i] );
						DrawFormatString(20,20*(sizeof(method)/sizeof(method[0])+1), GetColor( 0, 255, 0 ), "Now : %s", method[metn] );
					} else if(page == 3) {
						DrawFormatString( 0, 0, GetColor( 0, 255, 0 ), "Magnify : %s", buf );
						DrawFormatString( 0,20, GetColor( 0, 255, 0 ), "From 1 to 14000" );
						DrawFormatString(20,40, GetColor( 0, 255, 0 ), "Now : %d", mag );
					} else if(page == 4) {
						DrawFormatString( 0, 0, GetColor( 0, 255, 0 ), "Redraw(*dt[s]) : %s", buf );
						DrawFormatString( 0,20, GetColor( 0, 255, 0 ), "From 1 to 32767" );
						DrawFormatString(20,40, GetColor( 0, 255, 0 ), "Now : %d", draw );
						DrawFormatString(20,80, GetColor( 0, 255, 0 ), "dt[s] : %f", dt );
					} else if(page == 5) {
						DrawFormatString( 0, 0, GetColor( 0, 255, 0 ), "ε : %s", buf );
						if(eps<0.00000099 && eps!=0.0){
							DrawFormatString(20,20, GetColor( 0, 255, 0 ), "Now : %e", eps );
						} else {
							DrawFormatString(20,20, GetColor( 0, 255, 0 ), "Now : %f", eps );
						}
					}
					GetHitKeyStateAll( key );
				 	if(( key[KEY_INPUT_RETURN ] | key[KEY_INPUT_NUMPADENTER ] ) == 1 ) {
						if(page == 0) {
							viewtype = atoi(buf);
							center(viewtype);
						} else if(page == 1) {
							dt = atof(buf);
						} else if(page == 2) {
							metn = atoi(buf)%(sizeof(method)/sizeof(method[0]));
						} else if(page == 3) {
							mag = atoi(buf);
							if(mag > 14000) {
								mag = 14000;
							} else if(mag < 1) {
								mag = 1;
							}
						} else if(page == 4) {
							draw=atoi(buf);
							if(draw<=0)
								draw=1;
						} else if(page == 5) {
							eps = fabs(atof(buf));
						}
						break;
					}
					if( key[KEY_INPUT_X ] == 1) {
						break;
					}
					for(i=0;i<256;i++) {
						if(key[i] == 1) {
							tmpkey[i]++;
						} else {
							tmpkey[i] = 0;
						}
					}
					if( tmpkey[KEY_INPUT_LEFT ] == 1 )
					{
						page=(page+PAGE_MAX-1)%PAGE_MAX;
					}
					if( tmpkey[KEY_INPUT_RIGHT ] == 1 )
					{
						page=(page+1)%PAGE_MAX;
					}
					if( count < sizeof(buf)/sizeof(buf[0])-1 ) {
 						if( tmpkey[KEY_INPUT_NUMPAD0 ] == 1 || tmpkey[KEY_INPUT_0 ] == 1 ) {
							buf[count++] = '0';
							buf[count] = 0;
						}
						if( tmpkey[KEY_INPUT_NUMPAD1 ] == 1 || tmpkey[KEY_INPUT_1 ] == 1 ) {
							buf[count++] = '1';
							buf[count] = 0;
						}
						if( tmpkey[KEY_INPUT_NUMPAD2 ] == 1 || tmpkey[KEY_INPUT_2 ] == 1 ) {
							buf[count++] = '2';
							buf[count] = 0;
						}
						if( tmpkey[KEY_INPUT_NUMPAD3 ] == 1 || tmpkey[KEY_INPUT_3 ] == 1 ) {
							buf[count++] = '3';
							buf[count] = 0;
						}
						if( tmpkey[KEY_INPUT_NUMPAD4 ] == 1 || tmpkey[KEY_INPUT_4 ] == 1 ) {
							buf[count++] = '4';
							buf[count] = 0;
						}
						if( tmpkey[KEY_INPUT_NUMPAD5 ] == 1 || tmpkey[KEY_INPUT_5 ] == 1 ) {
							buf[count++] = '5';
							buf[count] = 0;
						}
						if( tmpkey[KEY_INPUT_NUMPAD6 ] == 1 || tmpkey[KEY_INPUT_6 ] == 1 ) {
							buf[count++] = '6';
							buf[count] = 0;
						}
						if( tmpkey[KEY_INPUT_NUMPAD7 ] == 1 || tmpkey[KEY_INPUT_7 ] == 1 ) {
							buf[count++] = '7';
							buf[count] = 0;
						}
						if( tmpkey[KEY_INPUT_NUMPAD8 ] == 1 || tmpkey[KEY_INPUT_8 ] == 1 ) {
							buf[count++] = '8';
							buf[count] = 0;
						}
						if( tmpkey[KEY_INPUT_NUMPAD9 ] == 1 || tmpkey[KEY_INPUT_9 ] == 1 ) {
							buf[count++] = '9';
							buf[count] = 0;
						}
						if( tmpkey[KEY_INPUT_SUBTRACT ] == 1 || tmpkey[KEY_INPUT_MINUS ] == 1 ) {
							buf[count++] = '-';
							buf[count] = 0;
						}
						if( tmpkey[KEY_INPUT_DECIMAL ] == 1 || tmpkey[KEY_INPUT_PERIOD ] == 1 ) {
							buf[count++] = '.';
							buf[count] = 0;
						} 
					}
					if( count > 0 && tmpkey[KEY_INPUT_BACK ] == 1 ) {
						buf[--count] = '\0';
					}
					// 裏画面の内容を表画面に反映
	       				ScreenFlip() ;
				}
				count = 0;
			}
			// 方向キーでカメラの座標を移動
			if( key[KEY_INPUT_LEFT ] == 1 )
			{
				CameraPos.x += (float)move[mn] * DirectX.x ;
				CameraPos.y += (float)move[mn] * DirectX.y ;
				CameraPos.z += (float)move[mn] * DirectX.z ;
				if((key[KEY_INPUT_LCONTROL] | key[KEY_INPUT_LCONTROL]) == 0) {
					ViewPoint.x += (float)move[mn] * DirectX.x ;
					ViewPoint.y += (float)move[mn] * DirectX.y ;
					ViewPoint.z += (float)move[mn] * DirectX.z ;
				}
			}
			if( key[KEY_INPUT_RIGHT ] == 1 )
			{
				CameraPos.x -= (float)move[mn] * DirectX.x ;
				CameraPos.y -= (float)move[mn] * DirectX.y ;
				CameraPos.z -= (float)move[mn] * DirectX.z ;
				if((key[KEY_INPUT_LCONTROL] | key[KEY_INPUT_LCONTROL]) == 0) {
					ViewPoint.x -= (float)move[mn] * DirectX.x ;
					ViewPoint.y -= (float)move[mn] * DirectX.y ;
					ViewPoint.z -= (float)move[mn] * DirectX.z ;
				}
			}
			if( key[KEY_INPUT_UP ] == 1 )
			{
				CameraPos.x += (float)move[mn] * DirectY.x ;
				CameraPos.y += (float)move[mn] * DirectY.y ;
				CameraPos.z += (float)move[mn] * DirectY.z ;
				if((key[KEY_INPUT_LCONTROL] | key[KEY_INPUT_LCONTROL]) == 0) {
					ViewPoint.x += (float)move[mn] * DirectY.x ;
					ViewPoint.y += (float)move[mn] * DirectY.y ;
					ViewPoint.z += (float)move[mn] * DirectY.z ;
				}
			}
			if( key[KEY_INPUT_DOWN ] == 1 )
			{
				CameraPos.x -= (float)move[mn] * DirectY.x ;
				CameraPos.y -= (float)move[mn] * DirectY.y ;
				CameraPos.z -= (float)move[mn] * DirectY.z ;
				if((key[KEY_INPUT_LCONTROL] | key[KEY_INPUT_LCONTROL]) == 0) {
					ViewPoint.x -= (float)move[mn] * DirectY.x ;
					ViewPoint.y -= (float)move[mn] * DirectY.y ;
					ViewPoint.z -= (float)move[mn] * DirectY.z ;
				}
			}
			if( key[KEY_INPUT_W ] == 1 )
			{
				CameraPos.x += (float)move[mn] * DirectZ.x ;
				CameraPos.y += (float)move[mn] * DirectZ.y ;
				CameraPos.z += (float)move[mn] * DirectZ.z ;
				if((key[KEY_INPUT_LCONTROL] | key[KEY_INPUT_LCONTROL]) == 0) {
					ViewPoint.x += (float)move[mn] * DirectZ.x ;
					ViewPoint.y += (float)move[mn] * DirectZ.y ;
					ViewPoint.z += (float)move[mn] * DirectZ.z ;
				}
			}
			if( key[KEY_INPUT_S ] == 1 )
			{
				CameraPos.x -= (float)move[mn] * DirectZ.x ;
				CameraPos.y -= (float)move[mn] * DirectZ.y ;
				CameraPos.z -= (float)move[mn] * DirectZ.z ;
				if((key[KEY_INPUT_LCONTROL] | key[KEY_INPUT_LCONTROL]) == 0) {
					ViewPoint.x -= (float)move[mn] * DirectZ.x ;
					ViewPoint.y -= (float)move[mn] * DirectZ.y ;
					ViewPoint.z -= (float)move[mn] * DirectZ.z ;
				}
			}
			for(i=0;i<256;i++) {
				if(key[i] == 1) {
					tmpkey[i]++;
				} else {
					tmpkey[i]=0;
				}
			}
			if( tmpkey[KEY_INPUT_Z ] == 1 || (tmpkey[KEY_INPUT_Z ] >= 60 && tmpkey[KEY_INPUT_Z ] % 60 == 0) )
			{
				mn=(mn+1)%(sizeof(move)/sizeof(move[0]));
			}
			if( tmpkey[KEY_INPUT_X ] == 1 || (tmpkey[KEY_INPUT_X ] >= 60 && tmpkey[KEY_INPUT_X ] % 30 == 0) )
			{
				mn=(mn-1+(sizeof(move)/sizeof(move[0])))%(sizeof(move)/sizeof(move[0]));
			}
			if( tmpkey[KEY_INPUT_R ] == 1 )
			{
				pn = 0;
				input();
				center(viewtype);
				wt = 0.0;
			}
			if( tmpkey[KEY_INPUT_H ] == 1 )
			{
				hide=(hide+2)%4;
			}
			if( tmpkey[KEY_INPUT_F ] == 1 )
			{
				while( ProcessMessage() == 0 && CheckHitKey( KEY_INPUT_ESCAPE ) == 0 ) {
					if(CheckHitKey(KEY_INPUT_G ) == 1)
						break;
					if(CheckHitKey(KEY_INPUT_R ) == 1) {
						pn = 0;
						input();
						center(viewtype);
						wt = 0.0;
						print_planets(x,wt,rot,move,CameraPos,rn,mn,hide);
						ScreenFlip() ;
					}
					
				}
			}
			if( tmpkey[KEY_INPUT_Q ] == 1 )
			{
				break;
			}
		}

	        // カメラの位置と注視点をセット
	        SetCameraPositionAndTargetAndUpVec( CameraPos, ViewPoint, UPVector ) ;
	
	        // 裏画面の内容を表画面に反映
	        ScreenFlip() ;
		c=(c+1)%draw;
	}

	DxLib_End() ;				// ＤＸライブラリ使用の終了処理

	return 0 ;				// ソフトの終了 
}
