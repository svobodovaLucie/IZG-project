/*!
 * @file
 * @brief This file contains implementation of gpu
 *
 * @author Tomáš Milet, imilet@fit.vutbr.cz
 */

#include <student/gpu.hpp>

#include <stdbool.h>
#include <stdio.h>
#include <iostream>
#include <cmath>

// Absolutni hodnota
#define ABS(x)          (((x) > 0) ? (x) : (-(x)))

// Minimum ze dvou zadanych hodnot
#define MIN(a, b)       (((a) < (b)) ? (a) : (b))

// Maximum ze dvou zadanych hodnot
#define MAX(a, b)       (((a) > (b)) ? (a) : (b))



#define t0 triangle.points[0]
#define t1 triangle.points[1]
#define t2 triangle.points[2]

#define tx0 triangle.points[0].gl_Position[0]
#define tx1 triangle.points[1].gl_Position[0]
#define tx2 triangle.points[2].gl_Position[0]

#define ty0 triangle.points[0].gl_Position[1]
#define ty1 triangle.points[1].gl_Position[1]
#define ty2 triangle.points[2].gl_Position[1]

#define tz0 triangle.points[0].gl_Position[2]
#define tz1 triangle.points[1].gl_Position[2]
#define tz2 triangle.points[2].gl_Position[2]

#define tw0 triangle.points[0].gl_Position[3]
#define tw1 triangle.points[1].gl_Position[3]
#define tw2 triangle.points[2].gl_Position[3]








struct Triangle{
  OutVertex points[3];
};


/*
struct Triangle {
  OutVertex point[3];
};

InVertex computeVertexID(GPUContext &ctx, InVertex inVertex, uint32_t i) {
    // kontrola, zda je zapnuté indexování
    if(ctx.vao.indexBuffer){  // vao.indexBuffer =! nullptr -> indexování je zapnuto
      // 32-bitové indexování
      if(ctx.vao.indexType == IndexType::UINT32)
        inVertex.gl_VertexID = ((uint32_t*)ctx.vao.indexBuffer)[i];
      // 16-bitové indexování
      else if(ctx.vao.indexType == IndexType::UINT16)
        inVertex.gl_VertexID = ((uint16_t*)ctx.vao.indexBuffer)[i];
      // 8-bitové indexování
      else if(ctx.vao.indexType == IndexType::UINT8 )
        inVertex.gl_VertexID = ((uint8_t *)ctx.vao.indexBuffer)[i];

    }else{  // vao.indexBuffer == nullptr -> indexování není zapnuto
      inVertex.gl_VertexID = i;
    }
    return inVertex;
}
*/
/*
InVertex readAtributes(GPUContext &ctx, InVertex inVertex) {
    for(uint32_t j = 0; j < maxAttributes; j++){
      auto const &currPos = ctx.vao.vertexAttrib[j];
      // pokud nejsou na dané pozici žádná data ke čtení, cyklus pokračuje dál
      if(!currPos.bufferData)
        continue;
      // pokud jsou na dané pozici nějaká data, zjistíme jaká a přečteme je z adresy buf_ptr + offset + stride*gl_VertexID
      if(currPos.type == AttributeType::FLOAT)
        inVertex.attributes[j].v1 = *(float*)(((uint8_t*)currPos.bufferData) + currPos.offset + currPos.stride * inVertex.gl_VertexID);
      if(currPos.type == AttributeType::VEC2 )
        inVertex.attributes[j].v2 = *(glm::vec2*)(((uint8_t*)currPos.bufferData) + currPos.offset + currPos.stride * inVertex.gl_VertexID);
      if(currPos.type == AttributeType::VEC3 )
        inVertex.attributes[j].v3 = *(glm::vec3*)(((uint8_t*)currPos.bufferData) + currPos.offset + currPos.stride * inVertex.gl_VertexID);
      if(currPos.type == AttributeType::VEC4 )
        inVertex.attributes[j].v4 = *(glm::vec4*)(((uint8_t*)currPos.bufferData) + currPos.offset + currPos.stride * inVertex.gl_VertexID);
    }
    return inVertex;
}
*/
/*
void runVertexAssembly(GPUContext &ctx, InVertex *inVertex, uint32_t v){
  computeVertexID(ctx, inVertex, v);
  readAtributes(ctx, inVertex);
}
*/


/*
void runPrimitiveAssembly(Triangle *triangle, GPUContext &ctx, uint32_t t) {
  for(uint32_t v = 0; v < 3; v++){
    InVertex inVertex;
    // triangle->points[0].gl_Position.x = 2;
    // runVertexAssembly(ctx, inVertex, t + v);

      inVertex = computeVertexID(ctx, inVertex, v);
      inVertex = readAtributes(ctx, inVertex);

        printf("before: x=%g, y=%g\n", triangle->point[0].gl_Position.x, triangle->point[0].gl_Position.y);

    ctx.prg.vertexShader(triangle->point[v], inVertex, ctx.prg.uniforms);
     printf("after: x=%g, y=%g\n", triangle->point[1].gl_Position.x, triangle->point[1].gl_Position.y);

 
  }
}
*/
/*
void createFragment(InFragment &inFragment, Triangle &triangle, glm::vec3 barycentric, glm::vec2 pixelCoord, GPUContext &ctx) {

}
*/
/*
float TriangleArea(float x0, float y0, float x1, float y1, float x2, float y2)
{
    float area_triangle;
    float a, b, c, s;

    a=std::sqrt(ABS((x0-x1)*(x0-x1)-(y0-y1)*(y0-y1)));
    b=std::sqrt(ABS((x1-x2)*(x1-x2)-(y1-y2)*(y1-y2)));
    c=std::sqrt(ABS((x0-x2)*(x0-x2)-(y0-y2)*(y0-y2)));

    s=(a+b+c)/2;

    area_triangle=std::sqrt(ABS((s*(s-a)*(s-b)*(s-c))));
    
    //printf("area: %g\n", area_triangle);

    return area_triangle;

}
*/
/*
void runPerspectiveDivision(Triangle *triangle) {
  for(unsigned i = 0; i < 3; i++) {
    triangle->points[i].gl_Position.x = triangle->points[i].gl_Position.x / triangle->points[i].gl_Position.w;
    triangle->points[i].gl_Position.y = triangle->points[i].gl_Position.y / triangle->points[i].gl_Position.w;
    triangle->points[i].gl_Position.z = triangle->points[i].gl_Position.z / triangle->points[i].gl_Position.w;
    // triangle->points[i].gl_Position.w = 1;
  }
}
*/
/*
void runViewportTransformation(Triangle *triangle, GPUContext &ctx) {
  for (unsigned i = 0; i < 3; i++) {
    triangle->points[i].gl_Position.x = (triangle->points[i].gl_Position.x * 0.5 + 0.5) * (ctx.frame.width);
    triangle->points[i].gl_Position.y = (triangle->points[i].gl_Position.y * 0.5 + 0.5) * (ctx.frame.height);
    // z a w zustava?
  }
}

*/

void runPerspectiveDivision(Triangle *triangle){
  for(int i = 0; i < 3; i++){
    triangle->points[i].gl_Position[0] /= triangle->points[i].gl_Position[3];
    triangle->points[i].gl_Position[1] /= triangle->points[i].gl_Position[3];
    triangle->points[i].gl_Position[2] /= triangle->points[i].gl_Position[3];
    //triangle->points[i].gl_Position[3] = 1;
  }
}

void runViewportTransformation(GPUContext &ctx, Triangle *triangle){
  for(int i = 0; i < 3; i++){
    triangle->points[i].gl_Position[0] += 1;
    triangle->points[i].gl_Position[1] += 1;
    triangle->points[i].gl_Position[0] *= (ctx.frame.width/2);
    triangle->points[i].gl_Position[1] *= (ctx.frame.height/2);
  }
}
/*
bool myBar(Triangle *triangle, u_int32_t x, u_int32_t y) {
  // vypocet obsahu

    //printf("before: x0=%g, y0=%g, x1=%g, y1=%g, x2=%g, y2=%g\n", triangle->points[0].gl_Position.x, triangle->points[0].gl_Position.y, triangle->points[1].gl_Position.x, triangle->points[1].gl_Position.y, triangle->points[2].gl_Position.x, triangle->points[2].gl_Position.y);
  // printf("UNOx: %d, DUEy: %d\n", x, y);
  float Sabc = TriangleArea(triangle->points[0].gl_Position.x, triangle->points[0].gl_Position.y, triangle->points[1].gl_Position.x, triangle->points[1].gl_Position.y, triangle->points[2].gl_Position.x, triangle->points[2].gl_Position.y);
  float Sabp = TriangleArea(triangle->points[0].gl_Position.x, triangle->points[0].gl_Position.y, triangle->points[1].gl_Position.x, triangle->points[1].gl_Position.y, x, y);
  float Sapc = TriangleArea(triangle->points[0].gl_Position.x, triangle->points[0].gl_Position.y, x, y, triangle->points[2].gl_Position.x, triangle->points[2].gl_Position.y);
  float Spbc = TriangleArea(x, y, triangle->points[1].gl_Position.x, triangle->points[1].gl_Position.y, triangle->points[2].gl_Position.x, triangle->points[2].gl_Position.y);

  float lambdaA = Spbc/Sabc;
  float lambdaB = Sapc/Sabc;
  float lambdaC = Sabp/Sabc;

  // printf("ABC: %f\n", Sabc);
  // printf("ABP: %f\n", Sabp);
  // printf("APC: %f\n", Sapc);
  // printf("PBC: %f\n", Spbc);
  // printf("LAMBDA: %f\n", lambdaA + lambdaB + lambdaC);

    // printf("before if : x=%g, y=%g\n", triangle->points[0].gl_Position.x, triangle->points[0].gl_Position.y);
    
  if (lambdaA + lambdaB + lambdaC == 1) {
    // bod nalezi trojuhelniku
    return true;
  }
  // bod nenalezi trojuhelniku
  return false;
}
*/
/*
void runBarycentric(Triangle *triangle, GPUContext &ctx) {
  // Nalezeni obalky (minX, maxX), (minY, maxY) trojuhleniku.
  int minX = MIN(MIN(triangle->points[0].gl_Position.x, triangle->points[1].gl_Position.x), triangle->points[2].gl_Position.x);
  int minY = MIN(MIN(triangle->points[0].gl_Position.y, triangle->points[1].gl_Position.y), triangle->points[2].gl_Position.y);
  int maxX = MAX(MAX(triangle->points[0].gl_Position.x, triangle->points[1].gl_Position.x), triangle->points[2].gl_Position.x);
  int maxY = MAX(MAX(triangle->points[0].gl_Position.y, triangle->points[1].gl_Position.y), triangle->points[2].gl_Position.y);

    // Oriznuti obalky (minX, maxX, minY, maxY) trojuhleniku podle rozmeru okna.
	minX = MAX(0, minX);
	minY = MAX(0, minY);
	maxX = MIN(ctx.frame.width - 1, maxX);
	maxY = MIN(ctx.frame.height - 1, maxY);

  // pruchod obalkou
	for (int y = minY; y <= maxY; y++) {
		bool even = (y - minY) % 2 == 0;

		int startX = even ? minX : maxX;
		int endX = even ? maxX : minX;
		int stepX = even ? 1 : -1;

		for (int x = startX; x != endX; x += stepX) {
      // cyklus pres vsechny hodnoty v obdelniku
      InFragment inFragment;
      OutFragment outFragment;
      if (myBar(triangle, x, y)) {
        ctx.prg.fragmentShader(outFragment, inFragment, ctx.prg.uniforms);
      }
		}
  }
}
/*
/*

void rasterizeTriangle(Triangle &triangle, GPUContext &ctx) {
  // spočítat hranice trojúhelníku

  // cyklus projde celý frame buffer
  //for(u_int32_t y = 0; y < ctx.frame.height; y++){
  //  for(u_int32_t x = 0; x < ctx.frame.width; x++){
      //if (ctx.frame.depth[(ctx.frame.width * y) + x]) {
          //InFragment inFragment;
          //runBarycentric(triangle, ctx);
          InFragment inFragment;
      OutFragment outFragment;
      //if (myBar(triangle, x, y)) {
        ctx.prg.fragmentShader(outFragment, inFragment, ctx.prg.uniforms);
          // createFragment(inFragment,primitive,barycentrics,pixelCoord,prg);
          //OutFragment outFragment;
          //ctx.prg.fragmentShader(outFragment,inFragment, ctx.prg.uniforms);
      //}
   // }
  //}
  
}
*/

InVertex computeVertexID(InVertex inVertex, GPUContext &ctx, int i){
  if(ctx.vao.indexBuffer != nullptr){
    if(ctx.vao.indexType == IndexType::UINT8){
      uint8_t *ind = (uint8_t*)ctx.vao.indexBuffer;
      inVertex.gl_VertexID = ind[i];
    }else if(ctx.vao.indexType == IndexType::UINT16){
      uint16_t *ind = (uint16_t*)ctx.vao.indexBuffer;
      inVertex.gl_VertexID = ind[i];
    }else{
      uint32_t *ind = (uint32_t*)ctx.vao.indexBuffer;
      inVertex.gl_VertexID = ind[i];
    }
  }else{
    inVertex.gl_VertexID = i;
  }
  return inVertex;
}

InVertex readAttributes(InVertex inVertex, GPUContext &ctx){
  for(int j = 0; j < maxAttributes; j++){
    if(ctx.vao.vertexAttrib[j].type != AttributeType::EMPTY){
      //Convert buffer to char*, so we can do pointer arithmetic by bytes
      char *buffer = (char*)ctx.vao.vertexAttrib[j].bufferData;
      int offset = ctx.vao.vertexAttrib[j].offset 
        + ctx.vao.vertexAttrib[j].stride * inVertex.gl_VertexID;

      if(ctx.vao.vertexAttrib[j].type == AttributeType::FLOAT){
        inVertex.attributes[j].v1 = *(float*)(buffer + offset);
      }
      else if(ctx.vao.vertexAttrib[j].type == AttributeType::VEC2){
        inVertex.attributes[j].v2 = *(glm::vec2*)(buffer + offset);
      }
      else if(ctx.vao.vertexAttrib[j].type == AttributeType::VEC3){
        inVertex.attributes[j].v3 = *(glm::vec3*)(buffer + offset);
      }
      else if(ctx.vao.vertexAttrib[j].type == AttributeType::VEC4){
        inVertex.attributes[j].v4 = *(glm::vec4*)(buffer + offset);
      }
    }
  }
  return inVertex;
}

void loadTriangle(GPUContext &ctx, uint32_t nofVertices, Triangle *triangle, int t){
  for(int i = t; i < t + 3; i++){//smyčka přes vrcholy
    InVertex inVertex; // vrchol, co vstupuje do vertex shader

    inVertex = computeVertexID(inVertex, ctx, i);
    inVertex = readAttributes(inVertex, ctx);

    printf("befor: t=%d: x=%g, y=%g\n", i, triangle->points[0].gl_Position.x, triangle->points[0].gl_Position.y);
    ctx.prg.vertexShader(triangle->points[i - t],inVertex,ctx.prg.uniforms);
    printf("after: t=%d: x=%g, y=%g\n", i, triangle->points[0].gl_Position.x, triangle->points[0].gl_Position.y);
  }
}

float dist(float p1[], float p2[]){
  return (float)sqrt(pow((double)(p2[0] - p1[0]), 2) + pow((double)(p2[1] - p1[1]), 2) * 1.0f);
}

float area(float a, float b, float c){
  float s = (a+b+c)/2;
  s = (s*(s-a)*(s-b)*(s-c));
  return (float)sqrt((double)((s < 0) ? -s : s));
}

glm::vec3 calcBary2(float p[], float a[], float b[], float c[]){
  float abc = area(dist(a, b), dist(b, c), dist(c, a));
  float bcp = area(dist(b, c), dist(c, p), dist(p, b));
  float acp = area(dist(a, c), dist(c, p), dist(p, a));
  float abp = area(dist(a, b), dist(b, p), dist(p, a));
  return glm::vec3(bcp/abc, acp/abc, abp/abc);
}

glm::vec3 getBary(Triangle *triangle, float x, float y){
  float p[4][2];
  for(int i = 0; i < 2; i++){
    p[0][i] = (float)triangle->points[0].gl_Position[i];
    p[1][i] = (float)triangle->points[1].gl_Position[i];
    p[2][i] = (float)triangle->points[2].gl_Position[i];
  }
  p[3][0] = x;
  p[3][1] = y;
  return glm::vec3(calcBary2(p[3], p[0], p[1], p[2]));
}

glm::vec3 interpolate3(float x, float y, Triangle triangle, int k){
    double area = tx0*(ty1-ty2) + tx1*(ty2-ty0) + tx2 * (ty0 - ty1);
    area = fabs((area/2.0));
    double area1 = x*(ty1-ty2) + tx1*(ty2-y) + tx2 * (y - ty1);
    area1=fabs(area1/2.0);

    double area2 = tx0*(y-ty2) + x*(ty2-ty0) + tx2 * (ty0 - y);
    area2=fabs(area2/2.0);

    double area3 = tx0*(ty1-y) + tx1*(y-ty0) + x * (ty0 - ty1);
    area3=fabs(area3/2.0);

    long double lambda1=area1/area;
    long double lambda2=area2/area;
    long double lambda3=area3/area;
    double s = lambda1/tw0 + lambda2/tw1 + lambda3/tw2;
    lambda1 = lambda1/(tw0*s);
    lambda2 = lambda2/(tw1*s);
    lambda3 = lambda3/(tw2*s);
    glm::vec3 q;
    double red = t0.attributes[k].v3[0]*lambda1 + t1.attributes[k].v3[0]*lambda2 + t2.attributes[k].v3[0]*lambda3;
    double green = t0.attributes[k].v3[1]*lambda1 + t1.attributes[k].v3[1]*lambda2 + t2.attributes[k].v3[1]*lambda3;
    double blue = t0.attributes[k].v3[2]*lambda1 + t1.attributes[k].v3[2]*lambda2 + t2.attributes[k].v3[2]*lambda3;
    q=glm::vec3(red,green,blue);
    return q;
}

glm::vec4 interpolate4(float x, float y, Triangle &triangle, int k){
    double area = tx0*(ty1-ty2) + tx1*(ty2-ty0) + tx2 * (ty0 - ty1);
    area = fabs((area/2.0));
    double area1 = x*(ty1-ty2) + tx1*(ty2-y) + tx2 * (y - ty1);
    area1=fabs(area1/2.0);

    double area2 = tx0*(y-ty2) + x*(ty2-ty0) + tx2 * (ty0 - y);
    area2=fabs(area2/2.0);

    double area3 = tx0*(ty1-y) + tx1*(y-ty0) + x * (ty0 - ty1);
    area3=fabs(area3/2.0);

    long double lambda1=area1/area;
    long double lambda2=area2/area;
    long double lambda3=area3/area;
    double s = lambda1/tw0 + lambda2/tw1 + lambda3/tw2;
    lambda1 = lambda1/(tw0*s);
    lambda2 = lambda2/(tw1*s);
    lambda3 = lambda3/(tw2*s);
    glm::vec4 q;
    double red = t0.attributes[k].v4[0]*lambda1 + t1.attributes[k].v4[0]*lambda2 + t2.attributes[k].v4[0]*lambda3;
    double green = t0.attributes[k].v4[1]*lambda1 + t1.attributes[k].v4[1]*lambda2 + t2.attributes[k].v4[1]*lambda3;
    double blue = t0.attributes[k].v4[2]*lambda1 + t1.attributes[k].v4[2]*lambda2 + t2.attributes[k].v4[2]*lambda3;
    double alpha = t0.attributes[k].v4[3]*lambda1 + t1.attributes[k].v4[3]*lambda2 + t2.attributes[k].v4[3]*lambda3;
    q=glm::vec4(red,green,blue,alpha);
    return q;
}

float is_in_triangle(float x, float y, Triangle &triangle, int number){
    
    double area = tx0*(ty1-ty2) + tx1*(ty2-ty0) + tx2 * (ty0 - ty1);
    area = fabs(area/2.0);
    double area1 = x*(ty1-ty2) + tx1*(ty2-y) + tx2 * (y - ty1);
    area1=(area1/2.0);

    double area2 = tx0*(y-ty2) + x*(ty2-ty0) + tx2 * (ty0 - y);
    area2=(area2/2.0);

    double area3 = tx0*(ty1-y) + tx1*(y-ty0) + x * (ty0 - ty1);
    area3=(area3/2.0);

    double total = area1+area2+area3;

    if (number==1){
        double lambda1=area1/area * tz0;
        double lambda2=area2/area * tz1;
        double lambda3=area3/area * tz2;

        return lambda1 + lambda2 + lambda3;

    }

    if(number==2){

    }
    if(area1>=0.0 && area2>= 0.0 && area3>= 0.0){
        return 1.0;
    }

    return 0;
}

void fixColor(OutFragment &outFragment){
    if (outFragment.gl_FragColor[0]<0){
        outFragment.gl_FragColor[0]=0.0f;
    }
    if (outFragment.gl_FragColor[1]<0){
        outFragment.gl_FragColor[1]=0.0f;
    }
    if (outFragment.gl_FragColor[2]<0){
        outFragment.gl_FragColor[2]=0.0f;
    }
    if (outFragment.gl_FragColor[0]>1){
        outFragment.gl_FragColor[0]=1.0f;
    }
    if (outFragment.gl_FragColor[1]>1){
        outFragment.gl_FragColor[1]=1.0f;
    }
    if (outFragment.gl_FragColor[2]>1){
        outFragment.gl_FragColor[2]=1.0f;
    }
}


void rasterize(GPUContext &ctx, Triangle *triangle) {
  // spocitat hranice trojuhlenika
    // Nalezeni obalky (minX, maxX), (minY, maxY) trojuhleniku.
  int minX = MIN(MIN(triangle->points[0].gl_Position.x, triangle->points[1].gl_Position.x), triangle->points[2].gl_Position.x);
  int minY = MIN(MIN(triangle->points[0].gl_Position.y, triangle->points[1].gl_Position.y), triangle->points[2].gl_Position.y);
  int maxX = MAX(MAX(triangle->points[0].gl_Position.x, triangle->points[1].gl_Position.x), triangle->points[2].gl_Position.x);
  int maxY = MAX(MAX(triangle->points[0].gl_Position.y, triangle->points[1].gl_Position.y), triangle->points[2].gl_Position.y);

    // Oriznuti obalky (minX, maxX, minY, maxY) trojuhleniku podle rozmeru okna.
	minX = MAX(0, minX);
	minY = MAX(0, minY);
	maxX = MIN(ctx.frame.width - 1, maxX);
	maxY = MIN(ctx.frame.height - 1, maxY);

  // pruchod obalkou
	for (int y = minY; y <= maxY; y++) {
		bool even = (y - minY) % 2 == 0;

		int startX = even ? minX : maxX;
		int endX = even ? maxX : minX;
		int stepX = even ? 1 : -1;

		for (int x = startX; x != endX; x += stepX) {

      // misko//////////////

                    if (is_in_triangle(x+0.5,y+0.5,*triangle,0)==1.0) {
                        InFragment inFragment;
                        inFragment.gl_FragCoord = glm::vec4((double)x + 0.5f, (double)y + 0.5f, is_in_triangle(x+0.5,y+0.5,*triangle,1), 1.0f);
                        for(size_t k=0; k< maxAttributes; k++){
                            if(ctx.prg.vs2fs[k]!=AttributeType::EMPTY){
                                if (ctx.prg.vs2fs[k]==AttributeType::FLOAT){

                                }
                                if (ctx.prg.vs2fs[k]==AttributeType::VEC2){

                                }
                                if (ctx.prg.vs2fs[k]==AttributeType::VEC3){
                                    inFragment.attributes[k].v3=interpolate3(x+0.5,y+0.5,*triangle,k);

                                }
                                if (ctx.prg.vs2fs[k]==AttributeType::VEC4){
                                    inFragment.attributes[k].v4=interpolate4(x+0.5,y+0.5,*triangle,k);
                                }
                            }
                        }


                        OutFragment outFragment;
                        ctx.prg.fragmentShader(outFragment, inFragment, ctx.prg.uniforms);
                        fixColor(outFragment);

                        if(inFragment.gl_FragCoord[2] < ctx.frame.depth[y*ctx.frame.width+x]){
                            *(ctx.frame.color+4*y*ctx.frame.width+x*4)   = ((float)*(ctx.frame.color+4*y*ctx.frame.width+x*4  )/255.f * (1.0f-outFragment.gl_FragColor[3])+outFragment.gl_FragColor[0]*outFragment.gl_FragColor[3])*255.f;
                            *(ctx.frame.color+4*y*ctx.frame.width+x*4+1) = ((float)*(ctx.frame.color+4*y*ctx.frame.width+x*4+1)/255.f * (1.0f-outFragment.gl_FragColor[3])+outFragment.gl_FragColor[1]*outFragment.gl_FragColor[3])*255.f;
                            *(ctx.frame.color+4*y*ctx.frame.width+x*4+2) = ((float)*(ctx.frame.color+4*y*ctx.frame.width+x*4+2)/255.f * (1.0f-outFragment.gl_FragColor[3])+outFragment.gl_FragColor[2]*outFragment.gl_FragColor[3])*255.f;
                            if (outFragment.gl_FragColor[3]>0.5f) {
                                *(ctx.frame.depth + y * ctx.frame.width + x) = inFragment.gl_FragCoord[2];
                            }

                        }
                    }

      // ////////////////
      // cyklus pres vsechny hodnoty v obdelniku
      /*InFragment inFragment;
      OutFragment outFragment;
      glm::vec3 bary = getBary(triangle, x+0.5, y+0.5);
      auto soucet_lambda = bary[0] + bary[1] + bary[2];
      // inFragment.gl_FragCoord = glm::vec4((double)x + 0.5f, (double)y + 0.5f, soucet_lambda, 1.0f);
       printf("Soucet lambda = %f\n", soucet_lambda);
      // if ((soucet_lambda >= 0.85 || soucet_lambda <= 1.15) || (bary[0] == 0 && (bary[1] + bary[2] == 1))) {
      //  if (soucet_lambda == 1.0) {
        ctx.prg.fragmentShader(outFragment, inFragment, ctx.prg.uniforms);
      //  }*/
		}
  }
}

//! [drawTrianglesImpl]
void drawTrianglesImpl(GPUContext &ctx,uint32_t nofVertices){
  /// \todo Tato funkce vykreslí trojúhelníky podle daného nastavení.<br>
  /// ctx obsahuje aktuální stav grafické karty.
  /// Parametr "nofVertices" obsahuje počet vrcholů, který by se měl vykreslit (3 pro jeden trojúhelník).<br>
  /// Bližší informace jsou uvedeny na hlavní stránce dokumentace.

  /*************** 1. Úkol - naprogramovat vertex assembly jednotku a pouštění vertex shaderu ***************/
  /*************** 2. Úkol - naprogramovat Primitive Assembly jednotku, rasterizaci a pouštění fragment shaderu ***************/

  // for(u_int32_t t = 0; t < nofVertices; t += 3) {
    // Triangle triangle;
    // runPrimitiveAssembly(&triangle, ctx, t);

Triangle triangle;
//printf("\n");
for(int i = 0; i < nofVertices; i += 3){
  loadTriangle(ctx, nofVertices, &triangle, i);
   // printf("t=%d: x=%g, y=%g\n", i, triangle.points[0].gl_Position.x, triangle.points[0].gl_Position.y);
   // printf("t=%d: x=%g, y=%g\n", i, triangle.points[1].gl_Position.x, triangle.points[1].gl_Position.y);
    //printf("t=%d: x=%g, y=%g\n", i, triangle.points[2].gl_Position.x, triangle.points[2].gl_Position.y);
  // printf("UNO: %f, DUE: %f\n", triangle->points[0].gl_Position.x, triangle->points[0].gl_Position.y);
  // printf("UNO: %f, DUE: %f\n", triangle->points[0].gl_Position.x, triangle->points[0].gl_Position.y);
  // printf("UNO: %f, DUE: %f\n", triangle->points[0].gl_Position.x, triangle->points[0].gl_Position.y);
  //printf("UNOx: %d, DUEy: %d\n", x, y);    
  
    runPerspectiveDivision(&triangle);
    runViewportTransformation(ctx, &triangle);
    //rasterizeTriangle(triangle, ctx);
    // runBarycentric(&triangle, ctx);
    rasterize(ctx, &triangle);

  }
}

//! [drawTrianglesImpl]

/**
 * @brief This function reads color from texture.
 *
 * @param texture texture
 * @param uv uv coordinates
 *
 * @return color 4 floats
 */
glm::vec4 read_texture(Texture const&texture,glm::vec2 uv){
  if(!texture.data)return glm::vec4(0.f);
  auto uv1 = glm::fract(uv);
  auto uv2 = uv1*glm::vec2(texture.width-1,texture.height-1)+0.5f;
  auto pix = glm::uvec2(uv2);
  //auto t   = glm::fract(uv2);
  glm::vec4 color = glm::vec4(0.f,0.f,0.f,1.f);
  for(uint32_t c=0;c<texture.channels;++c)
    color[c] = texture.data[(pix.y*texture.width+pix.x)*texture.channels+c]/255.f;
  return color;
}

/**
 * @brief This function clears framebuffer.
 *
 * @param ctx GPUContext
 * @param r red channel
 * @param g green channel
 * @param b blue channel
 * @param a alpha channel
 */
void clear(GPUContext&ctx,float r,float g,float b,float a){
  auto&frame = ctx.frame;
  auto const nofPixels = frame.width * frame.height;
  for(size_t i=0;i<nofPixels;++i){
    frame.depth[i] = 10e10f;
    frame.color[i*4+0] = static_cast<uint8_t>(glm::min(r*255.f,255.f));
    frame.color[i*4+1] = static_cast<uint8_t>(glm::min(g*255.f,255.f));
    frame.color[i*4+2] = static_cast<uint8_t>(glm::min(b*255.f,255.f));
    frame.color[i*4+3] = static_cast<uint8_t>(glm::min(a*255.f,255.f));
  }
}

