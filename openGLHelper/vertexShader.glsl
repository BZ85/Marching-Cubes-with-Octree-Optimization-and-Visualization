#version 150

in vec3 position;
in vec4 normal;




out vec4 col;

out vec3 viewPosition;
out vec3 viewNormal;


uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform mat4 normalMatrix;

void main()
{
  // compute the transformed and projected vertex position (into gl_Position) 
  // compute the vertex color (into col)
  //vec3 finalPosition;
 // gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0f);

  //col = color;
 

  col = normal;
 

 vec4 viewPosition4 = modelViewMatrix * vec4(position, 1.0f);
 viewPosition = viewPosition4.xyz;

 // final position in the normalized device coordinates space
 gl_Position = projectionMatrix * viewPosition4;

 // view-space normal
 viewNormal = normalize((normalMatrix * vec4(normal.xyz, 0.0f)).xyz);


}



