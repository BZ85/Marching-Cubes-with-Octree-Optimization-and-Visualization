#version 150


 in vec3 viewPosition;
 in vec3 viewNormal;

 in vec4 col;
 out vec4 c; // output color
 uniform vec4 La; // light ambient
 uniform vec4 Ld; // light diffuse
 uniform vec4 Ls; // light specular
 uniform vec3 viewLightDirection;
 uniform vec4 ka; // mesh ambient
 uniform vec4 kd; // mesh diffuse
 uniform vec4 ks; // mesh specular
 uniform float alpha; // shininess



void main()
{
 
// camera is at (0,0,0) after the modelview transformation

 vec3 eyedir = normalize(vec3(0, 0, 0) - viewPosition);
 
 vec3 reflectDir = -reflect(viewLightDirection, viewNormal);

 float d = max(dot(viewLightDirection, viewNormal), 0.0f);
 float s = max(dot(reflectDir, eyedir), 0.0f);

 c = ka * La + d * kd * Ld + pow(s, alpha) * ks * Ls;

 c = c * vec4(1, 0, 0, 1); // using the shading value to multiply the basic color

// c = vec4(viewNormal, 1);



}

