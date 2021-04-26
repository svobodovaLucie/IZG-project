/*!
 * @file
 * @brief This file contains implementation of gpu
 *
 * @author Tomáš Milet, imilet@fit.vutbr.cz
 */

#include <student/gpu.hpp>

void computeVertexID(GPUContext &ctx, InVertex &inVertex, uint32_t i) {
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
}

void readAtributes(GPUContext &ctx, InVertex &inVertex) {
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
}

void runVertexAssembly(GPUContext &ctx, InVertex &inVertex, uint32_t v){
  computeVertexID(ctx, inVertex, v);
  readAtributes(ctx, inVertex);
}

struct Triangle {
  OutVertex point[3];
};

void runPrimitiveAssembly(Triangle &triangle, GPUContext &ctx, uint32_t t) {

  // vrchol co vyleze z vertex shaderu (x, y, z, w)
  for(uint32_t v = 0; v < 3; v++){
    InVertex inVertex;
    runVertexAssembly(ctx, inVertex, t + v);
    ctx.prg.vertexShader(triangle.point[v], inVertex, ctx.prg.uniforms);
  }
}
/*
float area(glm::vec2 const&A,glm::vec2 const&B,glm::vec2 const&C){
  auto a=glm::length(B-A);
  auto b=glm::length(C-B);
  auto c=glm::length(A-C);
  auto s = (a+b+c)/2.f;
  return glm::sqrt(s*glm::abs(s-a)*glm::abs(s-b)*glm::abs(s-c));
}

glm::vec3 runPerspectiveDivision(Triangle &triangle, uint32_t width,uint32_t height,glm::vec2 x) {
  // clip-space
  glm::vec4 clip[3] = {
    triangle.point[0].gl_Position,
    triangle.point[1].gl_Position,
    triangle.point[2].gl_Position
  };

  // perspektivní dělení
  glm::vec3 ndc[3];
  for(int i=0;i<3;++i)
    ndc[i] = glm::vec3(clip[i])/clip[i].w;
  
  // viewport transformace
  glm::vec3 wc[3];
  for(int i=0;i<3;++i)
    wc[i] = glm::vec3((glm::vec2(ndc[i])*.5f+.5f)*glm::vec2(width, height),ndc[i].z);
  */
  /*
  float T = area(wc[0],wc[1],wc[2]);
  glm::vec3 l;
  for(int i=0;i<3;++i)l[i] = area(wc[(i+1)%3],wc[(i+2)%3],x)/T;
  return l;
  */
  //}

/*
  float T = area(wc[0],wc[1],wc[2]);
  glm::vec3 l;
  for(int i=0;i<3;++i)
    l[i] = area(wc[(i+1)%3],wc[(i+2)%3],x)/T;
  return l;
}

glm::vec3 barycentricsPerspective(Triangle &triangle,glm::vec2 x,uint32_t w,uint32_t h, glm::vec2 x){
  auto l = runPerspectiveDivision(&triangle,width,height, x);
  //glm::vec4 clip[3]={a.gl_Position,b.gl_Position,c.gl_Position};
  glm::vec4 clip[3] = {
    triangle.point[0].gl_Position,
    triangle.point[1].gl_Position,
    triangle.point[2].gl_Position
  };
  for(int i=0;i<3;++i)l[i] /= clip[i].w;
  l /= glm::dot(l,glm::vec3(1.f));
  return l;
}

void runViewportTransformation(Triangle &triangle, GPUContext &ctx) {

}

void rasterizeTriangle(Triangle &triangle, GPUContext &ctx) {

}
*/
//! [drawTrianglesImpl]
void drawTrianglesImpl(GPUContext &ctx,uint32_t nofVertices){
  //(void)ctx;
  //(void)nofVertices;
  /// \todo Tato funkce vykreslí trojúhelníky podle daného nastavení.<br>
  /// ctx obsahuje aktuální stav grafické karty.
  /// Parametr "nofVertices" obsahuje počet vrcholů, který by se měl vykreslit (3 pro jeden trojúhelník).<br>
  /// Bližší informace jsou uvedeny na hlavní stránce dokumentace.

  //OutVertex      outVertex;
  //InVertex  const inVertex;
  //Uniforms  const uniforms;


  /*************** 1. Úkol - naprogramovat vertex assembly jednotku a pouštění vertex shaderu ***************/
  VertexArray const &vao = ctx.vao;   // vertex array

  // for every vertex

  /*************** 2. Úkol - naprogramovat Primitive Assembly jednotku, rasterizaci a pouštění fragment shaderu ***************/
  
  for(u_int32_t t = 0; t < nofVertices; t += 3) {
    Triangle triangle;
    runPrimitiveAssembly(triangle, ctx, t);
    //runPerspectiveDivision(triangle, ctx.frame.width, ctx.frame.height);
  }
}
    //auto const&f:inFragments;
    //barycentricsPerspective(triangle, ctx.frame.width, ctx.frame.height, f.gl_FragCoord);
/*//////////////////////
float area(glm::vec2 const&A,glm::vec2 const&B,glm::vec2 const&C){
  auto a=glm::length(B-A);
  auto b=glm::length(C-B);
  auto c=glm::length(A-C);
  auto s = (a+b+c)/2.f;
  return glm::sqrt(s*glm::abs(s-a)*glm::abs(s-b)*glm::abs(s-c));
}

  glm::vec4 clip[3]={a.gl_Position,b.gl_Position,c.gl_Position};
  glm::vec3 ndc[3];
  for(int i=0;i<3;++i)ndc[i] = glm::vec3(clip[i])/clip[i].w;
  glm::vec3 wc[3];
  for(int i=0;i<3;++i)wc[i] = glm::vec3((glm::vec2(ndc[i])*.5f+.5f)*glm::vec2(w,h),ndc[i].z);
  float T = area(wc[0],wc[1],wc[2]);
  glm::vec3 l;
  for(int i=0;i<3;++i)l[i] = area(wc[(i+1)%3],wc[(i+2)%3],x)/T;

glm::vec3 barycentricsPerspective(OutVertex const&a,OutVertex const&b,OutVertex const&c,glm::vec2 x,uint32_t w,uint32_t h){
  auto l = barycentrics(a,b,c,x,w,h);
  glm::vec4 clip[3]={a.gl_Position,b.gl_Position,c.gl_Position};
  for(int i=0;i<3;++i)l[i] /= clip[i].w;
  l /= glm::dot(l,glm::vec3(1.f));
  return l;
}
*/ ////////////////////////////////
    //runViewportTransformation(triangle, ctx);


    //rasterizeTriangle(triangle, ctx);
    /*
    for(u_int32_t y = 0; y < ctx.frame.width; y++) {
      for(u_int32_t x = 0; x < ctx.frame.height; x++) {
        if () {
          
        }
      }
    }
    
  }
  */
 /************** GOOD
  for(u_int32_t t = 0; t < nofVertices; t++) {
    //runPrimitiveAssembly(triangle, ctx, t);
    InVertex inVertex;
    //runVertexAssembly(ctx, inVertex, t + v);
    //computeVertexID(ctx, inVertex, v);
    // kontrola, zda je zapnuté indexování
    if(ctx.vao.indexBuffer){  // vao.indexBuffer =! nullptr -> indexování je zapnuto
      // 32-bitové indexování
      if(ctx.vao.indexType == IndexType::UINT32)
        inVertex.gl_VertexID = ((uint32_t*)ctx.vao.indexBuffer)[t];
      // 16-bitové indexování
      else if(ctx.vao.indexType == IndexType::UINT16)
        inVertex.gl_VertexID = ((uint16_t*)ctx.vao.indexBuffer)[t];
      // 8-bitové indexování
      else if(ctx.vao.indexType == IndexType::UINT8 )
        inVertex.gl_VertexID = ((uint8_t *)ctx.vao.indexBuffer)[t];

    }else{  // vao.indexBuffer == nullptr -> indexování není zapnuto
      inVertex.gl_VertexID = t;
    }

    //readAtributes(ctx, inVertex);    
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
    OutVertex outVertex;
    ctx.prg.vertexShader(outVertex, inVertex, ctx.prg.uniforms);
  }
}
*/////////////////// GOOD END


    /**************/
    //InVertex inVertex;
    //runVertexAssembly(ctx, inVertex, t);
    //OutVertex outVertex;
    //ctx.prg.vertexShader(outVertex, inVertex, ctx.prg.uniforms);
    
    //runPerspectiveDivision(triangle, 100, 100);
    //runViewportTransformation(triangle, ctx);
    //rasterizeTriangle(triangle, ctx);
  
  
/****

void runPrimitiveAssembly(primitive,VertexArray vao,t,Program prg){
  for(every vertex v in triangle){
    InVertex inVertex;
    runVertexAssembly(inVertex,vao,t+v);
    prg.vertexShader(primitive.vertex,inVertex,prg.uniforms);
  }
}
 
void rasterizeTriangle(frame,primitive,prg){
  for(pixels in frame){
    if(pixels in primitive){
      InFragment inFragment;
      createFragment(inFragment,primitive,barycentrics,pixelCoord,prg);
      OutFragment outFragment;
      prg.fragmentShader(outFragment,inFragment,uniforms);
    }
  }
}
 
void drawTriangles(GPUContext&ctx,uint32_t nofVertices){
  for(every triangle t){
    Primitive primitive;
    runPrimitiveAssembly(primitive,ctx.vao,t,ctx.prg)
 
    runPerspectiveDivision(primitive)
    runViewportTransformation(primitive,ctx.frame)
    rasterizeTriangle(ctx.frame,primitive,ctx.prg);
  }
} 
 */



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

