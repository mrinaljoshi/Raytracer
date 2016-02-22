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

#define IMAGE_WIDTH		1500
#define IMAGE_HEIGHT	1500

int pixel_x= 0; int pixel_y=0;

Point E;	// = scene->camera->position;
Vec V;		// = scene->camera->viewDirection;
Flt f;// = scene->camera->focalDistance;
Vec U;// = scene->camera->orthoUp;
Flt phi;// = scene->camera->verticalFOV;*/
Point P; // point on the image corr to i,j pixel

typedef unsigned char u08;

SceneIO *scene = NULL;

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

/*
bool Shade(Point P_sphere, Vec N, L, ObjIO *intersected_obj) {
	// loop over all light source
	// find intersection with every object in this loop
	LightIO *var_light;
	var_light = scene->lights;
	 loop over all Light
		 loop over all object

		 if !intersect
			 save light and compute direct illu
		 else if intersect
			 check if distance > distance to source
				save light, compute 
				

	// Light sources 
	while (var_light->next != NULL)
	{
		Point P_light; // position of light source 
		Color Intensity;
		Intensity[0] = var_light->color[0];
		Intensity[1] = var_light->color[1];
		Intensity[2] = var_light->color[2];

		if (var_light->type == POINT_LIGHT)
		{
			Flt fatt;
			fatt = var_light->dropOffRate;
			P_light[0] = var_light->position[0];

			P_sphere
				N
				// Reflected ray = P_light - 
		}
		else if (var_light->type == DIRECTIONAL_LIGHT)
		{
			var_light->direction
		}

		var_light = scene->lights->next;
	}
}
*/
int main(int argc, char *argv[]) {
	Timer total_timer;
	total_timer.startTimer();

	loadScene("C:\\Users\\Mrinal\\Desktop\\Image Synthesis Assignment\\Assignment 2\\BasicRayTracer\\Scenes\\test2.scene");
	
	CImg<unsigned char> img(IMAGE_HEIGHT, IMAGE_WIDTH, 1, 3);        // Define a 640x400 color image with 8 bits per color component.
	                                
	for (int r = 0; r < IMAGE_HEIGHT; r++)
	{
		for (int c = 0; c < IMAGE_WIDTH; c++)
		{
			img(c, r, 0, 0) = 1;
			img(c, r, 0, 1) = 1;
			img(c, r, 0, 2) = 1;
		}
	}

	//img.display("My first CImg code");
	
	/*LightIO * light_source;
	light_source = scene->lights;
	Point P1;
		P1[0] = light_source->position[0];
		cout << P1[0];
	LightIO *save;
	save = light_source->next;
	Point P2;
	P2[0] = save->position[0];
	cout << P2[0];
	*/
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
			Flt *p;	//pointer to a float
			p = TraceRay(pixel_x, pixel_y, w, h, E, V, f, U, phi);
			P[0] = *(p); P[1] = *(p+1); P[2] = *(p+2);
		
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
					float t;

					if (IntersectSphere(P, C, R, E, t))
					{
						//call shadow ray (L. saved object pointer, point of interection on spehre, Normal)
						ObjIO *intersected_obj;
						intersected_obj = var_obj;

						MaterialIO *mat_obj;
						mat_obj = intersected_obj->material;

						Color Diff; Color Amb; Color Spec; Color Emissn; Flt Shin; Flt Trans;
						for (int i = 0; i < 2; i++) {
							Diff[i] = mat_obj->diffColor[i];
							Amb[i] = mat_obj->ambColor[i];
							Spec[i] = mat_obj->specColor[i];
							Emissn[i] = mat_obj->emissColor[i];
						}
						Shin = mat_obj->shininess;
						Trans = mat_obj->ktran;

						img(pixel_x, pixel_y, 0, 0) = 255 * (1 - Trans)*Amb[0] * Diff[0];
						img(pixel_x, pixel_y, 0, 1) = 255 * (1 - Trans)*Amb[0] * Diff[1];
						img(pixel_x, pixel_y, 0, 2) = 255 * (1 - Trans)*Amb[0] * Diff[2];

						Point P_sphere; //Point on the sphere where ray intersects
						Vec N; // Surface Normal at the point P_sphere
						for (int i = 0; i < 2; i++)
						{
							P_sphere[i] = E[i] + t*(P[i] - E[i]);
							N[i] = P_sphere[i] - C[i];
						}
						//bool Shade(Point P_sphere, Vec N , L, ObjIO *intersected_obj);
					}
				}
				/*
				else if(var_obj->type == POLYSET_OBJ)
				{
					PolySetIO *pset = (PolySetIO *)var_obj->data;
					//if (var_obj->type == POLYSET_TRI_MESH)
						PolygonIO *poly;
						VertexIO *vert;
						poly = pset->poly;

						for (int i = 0; i < pset->numPolys; i++, poly++) { //looping over all polygons
							vert = poly->vert;
							
							for (int j = 0; j < poly->numVertices; j++, vert++) {
								//Point Vertex0; Point Vertex1; Point Vertex2; 
								Flt polyvertices[9];
								
								for (int k = 0; k < 3; k++){
									polyvertices[3 * j + k] = vert->pos[k];

								}
								if (IntersectTriangle(polyvertices, P)) {
									ObjIO *intersected_obj;
									intersected_obj = var_obj;

									MaterialIO *mat_obj;
									mat_obj = intersected_obj->material;

									Color Diff; Color Amb; Color Spec; Color Emissn; Flt Shin; Flt Trans;
									for (int i = 0; i < 2; i++) {
										Diff[i] = mat_obj->diffColor[i];
										Amb[i] = mat_obj->ambColor[i];
										Spec[i] = mat_obj->specColor[i];
										Emissn[i] = mat_obj->emissColor[i];
									}
									Shin = mat_obj->shininess;
									Trans = mat_obj->ktran;

									img(pixel_x, pixel_y, 0, 0) = 255 * (1 - Trans)*Amb[0] * Diff[0];
									img(pixel_x, pixel_y, 0, 1) = 255 * (1 - Trans)*Amb[0] * Diff[1];
									img(pixel_x, pixel_y, 0, 2) = 255 * (1 - Trans)*Amb[0] * Diff[2];

								}
							}
						}
				}*/
				var_obj = var_obj->next;
				}
				}
			}
	img.display("some");
	/* write your ray tracer here */
	render();

	/* cleanup */
	if (scene != NULL) {
		deleteScene(scene);
	}

	total_timer.stopTimer();
	fprintf(stderr, "Total time: %.5lf secs\n\n", total_timer.getTime());
	
	return 1;
}
