#include <iostream>
#include "glm/glm.hpp"
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::mat3;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = SCREEN_WIDTH;
SDL_Surface* screen;
int t;
vector<Triangle> triangles;
vec3 cameraPos( 0, 0, -3.001 );
float f=SCREEN_WIDTH;
float yaw = 0.f;
mat3 R(cos(yaw),0,sin(yaw),0,1,0,-sin(yaw),0,cos(yaw));
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 lightPos(0,-0.5,-0.7);
vec3 lightPower = 14.0f*vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );
vec3 currentNormal;
vec3 currentReflectance;
float image_position_z = -4.f;

struct Pixel
{
    int x;
    int y;
    float zinv;
    vec3 pos3d;
};
struct Vertex
{
    vec3 position;
};
// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void Draw();
void VertexShader( const Vertex& v, Pixel& p );
void PixelShader( const Pixel& p );
//vec3 Illumination_old(const Vertex& v);
vec3 Illumination(const Pixel& v);
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );
void DrawLineSDL( SDL_Surface* surface, Pixel a, Pixel b, vec3 color );
void DrawPolygonEdges( const vector<vec3>& vertices );
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels,
                        vector<Pixel>& rightPixels );
void DrawPolygonRows( const vector<Pixel>& leftPixels,
              const vector<Pixel>& rightPixels );
void DrawPolygon( const vector<Vertex>& vertices );
bool InsideView(const int& x, const int& y);    //Function that determines whether the pixel is inside sdl                  surface or not.

int main( int argc, char* argv[] )
{
	LoadTestModel( triangles );
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	while( NoQuitMessageSDL() )
	{
        Update();
        Draw();
    }
	return 0;
}

void Update()
{
	Uint8* keystate = SDL_GetKeyState(0);

	if( keystate[SDLK_UP] ){
        //move camera forward:
        if (cameraPos.z<-2) {
            cameraPos.z+=0.05f;
            //Update focal length:
            f = (cameraPos.z - image_position_z)*SCREEN_HEIGHT;
        }
    }
		
	if( keystate[SDLK_DOWN] ){
        //move camera backward:
        if (cameraPos.z>image_position_z) {
            cameraPos.z-=0.05f;
            //Update focal length:
            f = (cameraPos.z - image_position_z)*SCREEN_HEIGHT;

        }
    }

	if( keystate[SDLK_RIGHT] ){
        if (yaw<M_PI/4) {
            yaw+=0.05;
            R = mat3(cos(yaw),0,sin(yaw),0,1,0,-sin(yaw),0,cos(yaw));
        }
    }

	if( keystate[SDLK_LEFT] ){
        if (yaw>-M_PI/4) {
            yaw-=0.05;
            R = mat3(cos(yaw),0,sin(yaw),0,1,0,-sin(yaw),0,cos(yaw));
        }
    }

	if( keystate[SDLK_w] ){
        lightPos.z+=0.05;
    }

	if( keystate[SDLK_s] ){
        lightPos.z-=0.05;
    }

	if( keystate[SDLK_d] ){
        lightPos.x+=0.05;
    }

	if( keystate[SDLK_a] ){
        lightPos.x-=0.05;
    }

	if( keystate[SDLK_e] ){
        lightPos.y+=0.05;
    }

	if( keystate[SDLK_q] ){
        lightPos.y-=0.05;
    }
}

vec3 Illumination(const Pixel& p){
    vec3 dir=glm::normalize(vec3(lightPos-p.pos3d));
    float radius=glm::length(p.pos3d-lightPos);
    float dot_prod = (glm::dot(dir, currentNormal));
    float max_val = (dot_prod > 0) ? dot_prod:0;
    
    float scale = max_val/(4*M_PI*pow(radius,2));
    glm::vec3 D =vec3(scale*lightPower);
    
    return currentReflectance*(D+indirectLightPowerPerArea);
}
 
bool InsideView(const int& x, const int& y){
    return (x>=0 && y>=0 && x<SCREEN_WIDTH && y<SCREEN_HEIGHT);
}

void PixelShader( const Pixel& p )
{
    if (InsideView(p.x, p.y)) {
        if(p.zinv>depthBuffer[p.y][p.x]){
            depthBuffer[p.y][p.x] = p.zinv;
            glm::vec3 color = Illumination(p);
            PutPixelSDL( screen, p.x, p.y, color );
        }
    }
}

void DrawPolygonRows( const vector<Pixel>& leftPixels,
                     const vector<Pixel>& rightPixels ){
    for (int i=0; i<leftPixels.size(); i++) {
        //Interpolate the depth over each row:
        size_t pixels=rightPixels[i].x-leftPixels[i].x+1;
       
        vector<Pixel> line( pixels );
        
        Interpolate( leftPixels[i], rightPixels[i], line );
        for (Pixel p: line) {
            PixelShader(p);
        } 
    }
}
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result ){
    int N = result.size();
    
    float step_x = (b.x-a.x) / float(max(N-1,1));
    float step_y = (b.y-a.y) / float(max(N-1,1));
    float step_zinv = (b.zinv-a.zinv) / float(max(N-1,1));
    vec3 step_3D = 1/float(max(N-1,1))*vec3((b.zinv*b.pos3d)-(a.zinv*a.pos3d));
    
    float current_x = a.x;
    float current_y = a.y;
    float current_zinv = a.zinv;
    
    vec3 current_3D = a.zinv*a.pos3d;
    Pixel current=a;
    
    for( int i=0; i<N; ++i )
    {
        result[i] = current;
        current_x += step_x;
        current_y += step_y;
        current_zinv += step_zinv;
        current_3D += step_3D;
        current.x = current_x + 0.5;
        current.y = current_y + 0.5;
        current.zinv = current_zinv;

        //multiply with z:
        current.pos3d = 1.f/current_zinv*current_3D;
    }
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels,
                        vector<Pixel>& rightPixels ){
    // 1. Find max and min y-value of the polygon
    // and compute the number of rows it occupies.
    int max_y = max(vertexPixels[0].y, vertexPixels[1].y);
    max_y = max(max_y, vertexPixels[2].y);
    
    int min_y = min(vertexPixels[0].y, vertexPixels[1].y);
    min_y=min(min_y, vertexPixels[2].y);
    size_t rows = max_y-min_y+1;
    
    // 2. Resize leftPixels and rightPixels
    leftPixels.resize(rows);
    rightPixels.resize(rows);
    
    //initialize:
    for( int i=0; i<rows; ++i )
    {
        leftPixels[i].x = +numeric_limits<int>::max();
        rightPixels[i].x = -numeric_limits<int>::max();
        leftPixels[i].y = rightPixels[i].y = i+min_y;
    }
    
    for( int i=0; i<vertexPixels.size(); ++i )
    {
        int j = (i+1)%vertexPixels.size();
        ivec2 delta = glm::abs(ivec2(vertexPixels[i].x,vertexPixels[i].y)-
                               ivec2(vertexPixels[j].x,vertexPixels[j].y));
        
        size_t pixels = glm::max( delta.x, delta.y ) + 1;
        vector<Pixel> line( pixels );
        
        Interpolate( vertexPixels[i], vertexPixels[j], line );
        
        for (Pixel elem: line) {
            if(leftPixels[elem.y-min_y].x>elem.x){
                leftPixels[elem.y-min_y].x=elem.x;
                leftPixels[elem.y-min_y].zinv=elem.zinv;
                leftPixels[elem.y-min_y].pos3d = elem.pos3d;
            }
            
            if (rightPixels[elem.y-min_y].x<elem.x) {
                rightPixels[elem.y-min_y].x=elem.x;
                rightPixels[elem.y-min_y].zinv=elem.zinv;
                rightPixels[elem.y-min_y].pos3d = elem.pos3d;
            }
        }
    }
}
void VertexShader(const Vertex& v, Pixel& p){
    //Step 1,
    //transform the vertices from the world coordinate system to the camera coordinate system.
    glm::vec3 P = v.position-cameraPos;
    P*=SCREEN_WIDTH;
    P= R*P;
    
    //Compute projection:
    int x = f*P.x/P.z;
    int y = f*P.y/P.z;
    
    //Translate to screen coordinate system:
    x+=SCREEN_WIDTH/2;
    y+=SCREEN_HEIGHT/2;
    
    float depth = sqrt(pow(cameraPos.x*SCREEN_WIDTH-P.x,2) +
                       pow(cameraPos.y*SCREEN_WIDTH-P.y,2) +
                       pow(cameraPos.z*SCREEN_WIDTH-P.z,2));
    //inverse of depth:
    depth = 1.0f/depth;
    
    p.x=x;
    p.y=y;
    p.zinv=depth;
    p.pos3d = v.position;
}

void DrawPolygon( const vector<Vertex>& vertices )
{
    int V = vertices.size();
    vector<Pixel> vertexPixels( V );
    
    for( int i=0; i<V; ++i ){
        VertexShader( vertices[i], vertexPixels[i] );
    }
    vector<Pixel> leftPixels(V);
    vector<Pixel> rightPixels(V);
    
    ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
    DrawPolygonRows( leftPixels, rightPixels );
}

void Draw()
{
	SDL_FillRect( screen, 0, 0 );

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);
	for( int y=0; y<SCREEN_HEIGHT; ++y )
        for( int x=0; x<SCREEN_WIDTH; ++x )
            depthBuffer[y][x] = 0.f;
	for( int i=0; i<triangles.size(); ++i )
	{
		vector<Vertex> vertices(3);
		vertices[0].position =triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;
        
        currentNormal = triangles[i].normal;
        currentReflectance = triangles[i].color;
        DrawPolygon(vertices);
	}
	
	if ( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}