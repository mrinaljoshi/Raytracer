#include <windows.h>
#include <stdio.h>
#include "scene_io.h"
#include "Timer.h"
#include <iostream> 
#include "Functions.h"
#include "Functions.cpp"
#include "CImg.h"

using namespace cimg_library;

using namespace std;

#define IMAGE_WIDTH		500
#define IMAGE_HEIGHT	500

int pixel_x= 0; int pixel_y=0;

Point E;	// = scene->camera->position;
Vec V;		// = scene->camera->viewDirection;
Flt f;// = scene->camera->focalDistance;
Vec U;// = scene->camera->orthoUp;
Flt phi;// = scene->camera->verticalFOV;*/
Point P; // point on the image corr to i,j pixel

typedef unsigned char u08;

SceneIO *scene = NULL;

typedef struct sortT {
Flt t;
Color Diff;
Color Amb; 
Color Spec;
Color Emissn;
Flt Shin; 
Flt Trans;
Point P_intersect;
Vec N_intersect;
};

static void loadScene(char *name) {
	/* load the scene into the SceneIO data structure using given parsing code */
	scene = readScene(name);

	/* hint: use the Visual Studio debugger ("watch" feature) to probe the
	   scene data structure and learn more about it for each of the given scenes */

	/* write any code to transfer from the scene data structure to your own here */
	/* */

	return;
}

/* just a place holder, feel free to edit */
void render(void) {
	int i, j, k;
	u08 *image = (u08 *)malloc(sizeof(u08) * IMAGE_HEIGHT * IMAGE_WIDTH * 3);
	u08 *ptr = image;

	for (j = 0; j < IMAGE_HEIGHT; j++) {
		for (i = 0; i < IMAGE_WIDTH; i++) {
			for (k = 0; k < 3; k++) {
				*(ptr++) = 0;
			}
		}
	}

	/* save out the image */
	/* */

	/* cleanup */
	free(image);

	return;
}

bool compare(sortT &first , sortT &second)
{
	return (first.t < second.t);
}

int main(int argc, char *argv[]) {
	Timer total_timer;
	total_timer.startTimer();

	loadScene("C:\\Users\\Mrinal\\Desktop\\Image Synthesis Assignment\\Assignment 2\\BasicRayTracer\\Scenes\\test4.scene");
	
	list <sortT> List_sort;
	
 
	CImg<unsigned char> img(IMAGE_HEIGHT, IMAGE_WIDTH, 1, 3);       // Define a 640x400 color image with 8 bits per color component.
	                                
	for (int r = 0; r < IMAGE_HEIGHT; r++)
	{
		for (int c = 0; c < IMAGE_WIDTH; c++)
		{
			img(c, r, 0, 0) = 1;
			img(c, r, 0, 1) = 1;
			img(c, r, 0, 2) = 1;
		}
	}
	
	for (int i = 0; i < 3; i++)
	{
		E[i] = scene->camera->position[i] ;
		V[i] = scene->camera->viewDirection[i];
		U[i] = scene->camera->orthoUp[i];
	}
	f = scene->camera->focalDistance;
	phi = scene->camera->verticalFOV;
	
	int w = IMAGE_WIDTH; int h = IMAGE_HEIGHT; // Image width & height 
	
	for (pixel_x = 0; pixel_x < w; pixel_x++) {
		for ( pixel_y = 0; pixel_y < h; pixel_y++) {
	//pixel_x = pixel_y = 300;
			Color Radiance;
			Radiance[0] = Radiance[1] = Radiance[2] = 0.0f;
			Flt *p;	//pointer to a float
			p = TraceRay(pixel_x, pixel_y, w, h, E, V, f, U, phi);
			P[0] = *(p); P[1] = *(p+1); P[2] = *(p+2); // P is point on the image plane corr to pixel_x, pixel_y
		
			ObjIO *var_obj;
			var_obj = scene->objects;
			
			while (var_obj != NULL)
			{
				
				if (var_obj->type == SPHERE_OBJ)
				{
					SphereIO *sphere = (SphereIO *)var_obj->data;

					Point C; //center of the sphere
					C[0] = sphere->origin[0];
					C[1] = sphere->origin[1];
					C[2] = sphere->origin[2];

					Flt R = sphere->radius; //radius of the sphere
					float t= 0.0f;

					if (IntersectSphere(P, C, R, E, &t))
					{
						ObjIO *intersected_obj;
						intersected_obj = var_obj;

						MaterialIO *mat_obj;
						mat_obj = intersected_obj->material;
						
						sortT comp;
						comp.t = t;
						
						//call shadow ray (L. saved object pointer, point of interection on spehre, Normal)
						//Color Diff; Color Amb; Color Spec; Color Emissn; Flt Shin; Flt Trans;
						for (int i = 0; i < 3; i++) {
							comp.Diff[i] = mat_obj->diffColor[i];
							comp.Amb[i] = mat_obj->ambColor[i];
							comp.Spec[i] = mat_obj->specColor[i];
							comp.Emissn[i] = mat_obj->emissColor[i];
						}
						comp.Shin = mat_obj->shininess;
						comp.Trans = mat_obj->ktran;
						Vec intermediate;
						Flt Nnorm = 0.0f; comp.P_intersect[0] = comp.P_intersect[1] = comp.P_intersect[2] = 0.0f;
						for (int i = 0; i < 3; i++)
						{
							comp.P_intersect[i] = E[i] + t*(P[i] - E[i]);
							
							//comp.N_intersect[i] =  (comp.P_intersect[i] - C[i])/(Nnorm);
						}
						Nnorm = sqrt((comp.P_intersect[0] - C[0])*(comp.P_intersect[0] - C[0]) + (comp.P_intersect[1] - C[1])*(comp.P_intersect[1] - C[1]) + (comp.P_intersect[2] - C[2])*(comp.P_intersect[2] - C[2]));
						for (int i = 0; i < 3; i++)
						{
							comp.N_intersect[i] =  (comp.P_intersect[i] - C[i])/(Nnorm);
						}
						List_sort.push_back(comp);

					}
				}
				
				else if(var_obj->type == POLYSET_OBJ)
				{
					PolySetIO *pset = (PolySetIO *)var_obj->data;
						PolygonIO *poly;
						VertexIO *vertex;
					
						poly = pset->poly;
						arr polyvertices;
						for (int i = 0; i < pset->numPolys; i++) { //looping over all polygons
							vertex = poly->vert;
							for (int j = 0; j < poly->numVertices; j++) {
								for (int k = 0; k < 3; k++) {
									polyvertices[3 * j + k] = vertex->pos[k];

								}
								vertex++;
							}
								Flt t = 0.0f;
								if (IntersectTriangle(polyvertices, P, E, &t)) {
									
									ObjIO *intersected_obj;
									intersected_obj = var_obj;

									MaterialIO *mat_obj;
									mat_obj = intersected_obj->material;
									
									sortT comp;
									comp.t = t;

									//Color Diff; Color Amb; Color Spec; Color Emissn; Flt Shin; Flt Trans;
									for (int i = 0; i < 3; i++) {
										comp.Diff[i] = mat_obj->diffColor[i];
										comp.Amb[i] = mat_obj->ambColor[i];
										comp.Spec[i] = mat_obj->specColor[i];
										comp.Emissn[i] = mat_obj->emissColor[i];
									}
									comp.Shin = mat_obj->shininess;
									comp.Trans = mat_obj->ktran;

									for (int i = 0; i < 3; i++)
									{
										comp.P_intersect[i] = E[i] + t*(P[i] - E[i]);
										//comp.N_intersect[i] = comp.P_intersect[i] - C[i];
									}

									Vec v0; Vec v1; Vec v2;
									for (int i = 0; i < 3; i++)
									{
										v0[i] = polyvertices[i];
										v1[i] = polyvertices[i + 3];
										v2[i] = polyvertices[i + 6];
									}
									Vec v10; Vec v20; Vec N;
									for (int i = 0; i < 3; i++)
									{
										v10[i] = v1[i] - v0[i];
										v20[i] = v2[i] - v0[i];
									}
									Point N_tri;
									N_tri[0] = v10[1]*v20[2] - v10[2] * v20[1];
									N_tri[1] = v10[2]*v20[0] - v10[0] * v20[2];
									N_tri[2] = v10[0]*v20[1] - v10[1] * v20[0];
									
									Flt Nnorm;
									Nnorm= sqrt(N_tri[0] * N_tri[0] + N_tri[1] * N_tri[1] + N_tri[2] * N_tri[2]);
									comp.N_intersect[0] = N_tri[0]/Nnorm;
									comp.N_intersect[1] = N_tri[1] / Nnorm;
									comp.N_intersect[2] = N_tri[2] / Nnorm;

									List_sort.push_back(comp);
								}
							poly++;
						}
				}
				var_obj = var_obj->next;
				}
					
				
				if (!List_sort.empty())
					{
						List_sort.sort(compare); // found the nearest t. 
						sortT k = *List_sort.begin();

						Point P_shadow; // Point of intersection plus epsilon
						for (int i = 0; i < 3; i++) {
							P_shadow[i] = k.P_intersect[i] + 0.001*(k.N_intersect[i]);
						}
						
						LightIO *light_source;
						light_source = scene->lights;
						Point P_light; Vec RayDir; Flt fatt; Flt tmax;// dist from P_shadow to light source
					
						while (light_source != NULL)
						{
							Vec S; Vec S_total;
							S[0] =S[1]=S[2] = 1.0f;
							S_total[0] = S_total[1] = S_total[2] = 0.0f;
							if (light_source->type == POINT_LIGHT){//light_source->type == POINT_LIGHT) {
								
								for (int i = 0; i < 3; i++) {
									P_light[i] = light_source->position[i];
									RayDir[i] = P_light[i] - P_shadow[i];
									}
								tmax = sqrt((P_shadow[0] - P_light[0])*(P_shadow[0] - P_light[0]) + ((P_shadow[1] - P_light[1])*(P_shadow[1] - P_light[1])) + ((P_shadow[2] - P_light[2])*(P_shadow[2] - P_light[2])));

								fatt = min(1.0f, (1.0f/(0.25f + 0.1f*tmax + 0.01f*tmax*tmax)) );
								}
							else if (light_source->type == DIRECTIONAL_LIGHT) {
								for (int i = 0; i < 3; i++) {
									RayDir[i] = light_source->direction[i];//light_source->direction[i];
									}
								fatt = 1.0f;
								//tmax = INFINITE;
								//tmax = 1000000000000.0f;
								}

							Flt RayDirNorm = sqrt(RayDir[0] * RayDir[0] + RayDir[1] * RayDir[1] + RayDir[2] * RayDir[2]);
							
							// compute parameters for that light source and object of intersection.
							Vec I; // illumination due to light source
							I[0] = light_source->color[0]; I[1] = light_source->color[1]; I[2] = light_source->color[2];
							
							//I[0] = scene->lights->color[0]; I[1] = scene->lights->color[1]; I[2] = scene->lights->color[2];

							Vec Refl; // reflection vector
							Flt NdotL;
							
							RayDir[0] = RayDir[0] / RayDirNorm; RayDir[1] = RayDir[1] / RayDirNorm; RayDir[2] = RayDir[2] / RayDirNorm;

							NdotL = (k.N_intersect[0] * RayDir[0] + k.N_intersect[1] * RayDir[1] + k.N_intersect[2] * RayDir[2]);
							for (int i = 0; i < 3; i++) {
								Refl[i] = 2 * NdotL*k.N_intersect[i] - RayDir[i];
							}
							Flt RdotV; Vec View;
							Flt Viewnorm;
							View[0] = E[0] - P_shadow[0]; View[1] = E[1] - P_shadow[1]; View[2] = E[2] - P_shadow[2];

							Viewnorm = sqrt(View[0] * View[0] + View[1] * View[1] + View[2] * View[2]);

							RdotV = (Refl[0] * View[0] + Refl[1] * View[1] + Refl[2] * View[2])/Viewnorm;

							    ObjIO *shadow_obj;
								shadow_obj = scene->objects;
									while (shadow_obj != NULL) {

									if (shadow_obj->type == SPHERE_OBJ)
									{
										SphereIO *sphere = (SphereIO *)shadow_obj->data;

										Point C; //center of the sphere
										C[0] = sphere->origin[0];
										C[1] = sphere->origin[1];
										C[2] = sphere->origin[2];

										Flt R = sphere->radius; //radius of the sphere
										float t = 0.0f;
										//Vec RayDir;
										
										if (IntersectSphere_shadow(RayDir, C, R, P_shadow, &t)) {
											ObjIO *intersect_obj;
											intersect_obj = shadow_obj;
											Flt Cdmax; 
											Cdmax = max(intersect_obj->material->diffColor[0], max(intersect_obj->material->diffColor[1], intersect_obj->material->diffColor[2]));
											//cout << intersect_obj->material->diffColor[0] << endl << "diff0"<< endl<<intersect_obj->material->diffColor[1]<<endl<<"diff1"<< intersect_obj->material->diffColor[2]<<endl<<"diff2"<<endl<<Cdmax<<endl;
											if (light_source->type == POINT_LIGHT)
											{
												if (t < tmax && t>0.0001) {
													if (intersect_obj->material->ktran == 0)
													{
														S[0] = S[1] = S[2] = 0.0f;
														break;
													}

													else {

														for (int i = 0; i < 3; i++) {
															S[i] = (intersect_obj->material->ktran*intersect_obj->material->diffColor[i] * S[i]) / Cdmax;
														}
													}
												}
											}
											else if (light_source->type == DIRECTIONAL_LIGHT)
											{
												if (t>0.0001) {
													if (intersect_obj->material->ktran == 0)
													{
														S[0] = S[1] = S[2] = 0.0f;
														break;
													}

													else {

														for (int i = 0; i < 3; i++) {
															S[i] = (intersect_obj->material->ktran*intersect_obj->material->diffColor[i] * S[i]) / Cdmax;
														}
													}
												}
											}
											
										}
										}
									
									else if (shadow_obj->type == POLYSET_OBJ)
									{
										PolySetIO *pset = (PolySetIO *)shadow_obj->data;
										PolygonIO *poly;
										VertexIO *vertex;

										poly = pset->poly;
										arr polyvertices;
										for (int i = 0; i < pset->numPolys; i++) { //looping over all polygons
											vertex = poly->vert;
											for (int j = 0; j < poly->numVertices; j++) {
												for (int k = 0; k < 3; k++) {
													polyvertices[3 * j + k] = vertex->pos[k];

												}
												vertex++;
											}
											Flt t = 0.0f;
											if (IntersectTriangle_shadow(polyvertices, P_shadow, RayDir, &t)) {
												ObjIO *intersect_obj;
												intersect_obj = shadow_obj;
												Flt Cdmax; Cdmax = max(intersect_obj->material->diffColor[0], max(intersect_obj->material->diffColor[1], intersect_obj->material->diffColor[2]));;
												
												if (light_source->type == POINT_LIGHT)
												{
													if (t < tmax && t>0.0001) {
														if (intersect_obj->material->ktran == 0)
														{
															S[0] = S[1] = S[2] = 0.0f;
															break;
														}

														else {

															for (int i = 0; i < 3; i++) {
																S[i] = (intersect_obj->material->ktran*intersect_obj->material->diffColor[i] * S[i]) / Cdmax;
															}
														}
													}
												}
												else if (light_source->type == DIRECTIONAL_LIGHT)
												{
													if (t>0.0001) {
														if (intersect_obj->material->ktran == 0)
														{
															S[0] = S[1] = S[2] = 0.0f;
															break;
														}

														else {

															for (int i = 0; i < 3; i++) {
																S[i] = (intersect_obj->material->ktran*intersect_obj->material->diffColor[i] * S[i]) / Cdmax;
															}
														}
													}
												}
												
												
											}
											poly++;
										}
									}
									
									shadow_obj = shadow_obj->next;
								}
									
									Flt shiny;
									
										if (RdotV > 0 && k.Shin > 0) { 
										shiny = pow(RdotV, k.Shin);
									}
									else { 
										shiny = 0.0f; 
									}
									//cout << "pow" << shiny << endl;
							for (int i = 0; i < 3; i++) {
								Radiance[i] = Radiance[i] + S[i] * fatt * I[i] * ((1.0f - k.Trans)*k.Diff[i] * max(NdotL, 0.0f) + k.Spec[i] * shiny);
							}

							light_source = light_source->next;
						}

						//  * fatt*I[i] * ((1.0f - k.Trans)*k.Diff[i] * max(NdotL, 0) + k.Spec[i] * pow(RdotV, k.Shin)))
						img(pixel_x, IMAGE_HEIGHT - pixel_y, 0, 0) = 255 * min(1.0f,(Radiance[0]+(1.0f - k.Trans)*k.Amb[0] * k.Diff[0]));
						img(pixel_x, IMAGE_HEIGHT - pixel_y, 0, 1) = 255 * min(1.0f,(Radiance[1]+(1.0f - k.Trans)*k.Amb[0] * k.Diff[1]));
						img(pixel_x, IMAGE_HEIGHT - pixel_y, 0, 2) = 255 * min(1.0f,(Radiance[2] +(1.0f - k.Trans)*k.Amb[0] * k.Diff[2]));
					}
					List_sort.clear();
				}  // pixel_y
			} // pixel_x
	img.display("Image");
	/* write your ray tracer here */
	render();

	/* cleanup */
	if (scene != NULL) {
		deleteScene(scene);
	}

	total_timer.stopTimer();
	//fprintf(stderr, "Total time: %.5lf secs\n\n", total_timer.getTime());
	
	return 1;
}	
