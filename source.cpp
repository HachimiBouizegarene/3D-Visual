

#include <fstream>
#include <strstream>
#include <algorithm>
#include "olcConsoleGameEngineGL.h"
using namespace std;


struct vec3d
{
	float x, y, z;
	float w = 1.0f;
};

struct vec2d {
	float x, y;
};

struct triangle
{
	vec3d p[3];
	vec2d t[3];
	CHAR_INFO infos;
};

struct mesh
{
	vector<triangle> tris;

	bool LoadFromObjectFile(string sFilename)
	{
		ifstream f(sFilename);
		if (!f.is_open())
			return false;


		vector<vec3d> verts;

		while (!f.eof())
		{
			char line[128];
			f.getline(line, 128);

			strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				vec3d v;
				s >> junk >> v.x >> v.y >> v.z;
				verts.push_back(v);
			}

			if (line[0] == 'f')
			{
				int f[3];
				s >> junk >> f[0] >> f[1] >> f[2];
				tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
			}
		}

		return true;
	}
};

struct mat4x4
{
	float m[4][4] = { 0 };
};

class olcEngine3D : public olcConsoleGameEngine
{
public:
	olcEngine3D()
	{
		m_sAppName = L"3D-Visual";
	}


private:
	mesh meshCube;
	mat4x4 matProj;
	mat4x4 matRot;
	float fTheta;
	olcSprite* sprTex1;
	vec3d vCamera;
	vec3d lookCamera;
	vec3d		lookCameraRotated;
	vec3d lookCameraRight;
	vec3d lookCameraRightRotated;

	float fYaw;
public:

	vec3d Vector_IntersectPlane(vec3d& plane_p, vec3d& plane_n, vec3d& lineStart, vec3d& lineEnd,float &t)
	{

		vec3d vector_director = Vector_sub(lineStart, lineEnd); 
		float a = Vector_dot(Vector_sub(plane_p, lineStart), plane_n);
		float b = Vector_dot(vector_director, plane_n);
		t = a / b;
		vec3d ret = Vector_add(lineStart, Vector_mul(vector_director, t));
		return ret;
	}

	int Triangle_ClipAgainstPlane(vec3d plane_p, vec3d plane_n, triangle& in_tri, triangle& out_tri1, triangle& out_tri2)
	{

		plane_n = Vector_normalize(plane_n);

	
		auto dist = [&](vec3d& p)
		{
			vec3d n = Vector_normalize(p);
			return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_dot(plane_n, plane_p));
		};


		vec3d* inside_points[3];  int nInsidePointCount = 0;
		vec3d* outside_points[3]; int nOutsidePointCount = 0;
		vec2d* inside_tex[3]; int nInsideTexCount = 0;
		vec2d* outside_tex[3]; int nOutsideTexCount = 0;

		float d0 = dist(in_tri.p[0]);
		float d1 = dist(in_tri.p[1]);
		float d2 = dist(in_tri.p[2]);

		if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; inside_tex[nInsideTexCount++] = &in_tri.t[0]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; outside_tex[nOutsideTexCount++] = &in_tri.t[0]; }
		if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; inside_tex[nInsideTexCount++] = &in_tri.t[1];
		}
		else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; outside_tex[nOutsideTexCount++] = &in_tri.t[1];
		}
		if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; inside_tex[nInsideTexCount++] = &in_tri.t[2];
		}
		else { outside_points[nOutsidePointCount++] = &in_tri.p[2];outside_tex[nOutsideTexCount++] = &in_tri.t[2];
		}



		if (nInsidePointCount == 0)
		{


			return 0; 
		}

		if (nInsidePointCount == 3)
		{

			out_tri1 = in_tri;

			return 1; 
		}

		if (nInsidePointCount == 1 && nOutsidePointCount == 2)
		{

			out_tri1.infos.Attributes = in_tri.infos.Attributes;
			out_tri1.infos.Char = in_tri.infos.Char;


			out_tri1.p[0] = *inside_points[0];
			out_tri1.t[0] = *inside_tex[0];


			float t;
			out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
			out_tri1.t[1].x = t * (outside_tex[0]->x - inside_tex[0]->x) + inside_tex[0]->x;
			out_tri1.t[1].y = t * (outside_tex[0]->y - inside_tex[0]->y) + inside_tex[0]->y;

			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1], t);

			out_tri1.t[2].x = t * (outside_tex[1]->x - inside_tex[0]->x) + inside_tex[0]->x;
			out_tri1.t[2].y = t * (outside_tex[1]->y - inside_tex[0]->y) + inside_tex[0]->y;

			return 1; 
		}

		if (nInsidePointCount == 2 && nOutsidePointCount == 1)
		{

			out_tri1.infos.Attributes = in_tri.infos.Attributes;
			out_tri1.infos.Char = in_tri.infos.Char;

			out_tri2.infos.Attributes = in_tri.infos.Attributes;
			out_tri2.infos.Char = in_tri.infos.Char;


			out_tri1.p[0] = *inside_points[0];
			out_tri1.p[1] = *inside_points[1];
			out_tri1.t[0] = *inside_tex[0];
			out_tri1.t[1] = *inside_tex[1];

			float t;
			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0],t);

			out_tri1.t[2].x = t * (outside_tex[0]->x - inside_tex[0]->x) + inside_tex[0]->x;
			out_tri1.t[2].y = t * (outside_tex[0]->y - inside_tex[0]->y) + inside_tex[0]->y;

			out_tri2.p[0] = *inside_points[1];
			out_tri2.p[1] = out_tri1.p[2];
			out_tri2.t[0] = *inside_tex[1];
			out_tri2.t[1] = out_tri1.t[1];
			out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0],t);
			out_tri2.t[2].x = t * (outside_tex[0]->x - inside_tex[1]->x) + inside_tex[1]->x;
			out_tri2.t[2].y = t * (outside_tex[0]->y - inside_tex[1]->y) + inside_tex[1]->y;
			return 2; 
		}
	}

	mat4x4 Matrix_MakeRotationZ(float fAngleRad)
	{
		mat4x4 matrix;
		matrix.m[0][0] = cosf(fAngleRad);
		matrix.m[0][1] = sinf(fAngleRad);
		matrix.m[1][0] = -sinf(fAngleRad);
		matrix.m[1][1] = cosf(fAngleRad);
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		return matrix;
	}
	mat4x4 Matrix_MakeRotationY(float fAngleRad)
	{
		mat4x4 matrix;
		matrix.m[0][0] = cosf(fAngleRad);
		matrix.m[0][2] = sinf(fAngleRad);
		matrix.m[2][0] = -sinf(fAngleRad);
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = cosf(fAngleRad);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}

	mat4x4 Matrix_MakeRotationX(float fAngleRad)
	{
		mat4x4 matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = cosf(fAngleRad);
		matrix.m[1][2] = sinf(fAngleRad);
		matrix.m[2][1] = -sinf(fAngleRad);
		matrix.m[2][2] = cosf(fAngleRad);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}

	mat4x4 Matrix_MakeIdentity() {
		mat4x4 matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = 1.0f;
		matrix.m[2][3] = 1.0f;

		return matrix;
	}
	vec3d MultiplyMatrixVector(vec3d& v, mat4x4 mat) {
		vec3d ret;


		ret.x = v.x * mat.m[0][0] + v.y * mat.m[1][0] + v.z * mat.m[2][0] + v.w * mat.m[3][0];
		ret.y = v.x * mat.m[0][1] + v.y * mat.m[1][1] + v.z * mat.m[2][1] + v.w * mat.m[3][1];
		ret.z = v.x * mat.m[0][2] + v.y * mat.m[1][2] + v.z * mat.m[2][2] + v.w * mat.m[3][2];
		ret.w = v.x * mat.m[0][3] + v.y * mat.m[1][3] + v.z * mat.m[2][3] + v.w * mat.m[3][3];

		return ret;
	}

	vec3d Vector_div(vec3d v, float d) {
		return {
		v.x / d,
		v.y / d,
		v.z / d,
		v.w
		};
	}

	vec3d Vector_mul(vec3d v, float d) {
		return {
		v.x * d,
		v.y * d,
		v.z * d,
		v.w
		};
	}

	CHAR_INFO GetColour(float lum)
	{



		short bg_col, fg_col;
		wchar_t sym;
		int pixel_bw = (int)(12.0f * lum);
		switch (pixel_bw)
		{
		case 0: bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID; break;

		case 1: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_QUARTER; break;
		case 2: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_HALF; break;
		case 3: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 4: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_SOLID; break;

		case 5: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_QUARTER; break;
		case 6: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_HALF; break;
		case 7: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 8: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_SOLID; break;

		case 9:  bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_QUARTER; break;
		case 10: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_HALF; break;
		case 11: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_THREEQUARTERS; break;
		case 12: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_SOLID; break;
		default:
			bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID;
		}

		CHAR_INFO c;
		c.Attributes = bg_col | fg_col;
		c.Char.UnicodeChar = sym;
		return c;
	}

	float Vector_dot(vec3d v1, vec3d v2) {
		float ret = v1.x * v2.x;
		ret += v1.y * v2.y;
		ret += v1.z * v2.z;
		return ret;
	}

	vec3d Vector_toLocal(vec3d v, vec3d pos , vec3d forward) {

		//a corriger pour permettre une rotation autour de X et Z
		vec3d up = { 0.0f, 1.0f, 0.0f, 1.0f };
		vec3d a = Vector_mul(forward, Vector_dot(up, forward));

		vec3d newUp = Vector_sub(up, a);

		newUp = Vector_normalize(newUp);
		vec3d right = Vector_crossProduct(newUp, forward);

		//vecteur qui va du point v a la position de la camera 
		vec3d vecPos = Vector_sub(v, pos);

		vec3d ret;
		ret.x = Vector_dot(vecPos, right);
		ret.y= Vector_dot(vecPos, newUp);
		ret.z = Vector_dot(vecPos, forward);

		return ret;
	}

	vec3d Vector_crossProduct(vec3d v1, vec3d v2) {
		vec3d ret;

		ret.x = v1.y * v2.z - v1.z * v2.y;
		ret.y = v1.z * v2.x - v1.x * v2.z;
		ret.z = v1.x * v2.y - v1.y * v2.x;

		return ret;
	}

	vec3d Vector_sub(vec3d v1, vec3d v2) {
		return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w };
	}

	vec3d Vector_add(vec3d v1, vec3d v2) {
		return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w };
	}
	float Vector_length(vec3d v) {
		return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	}

	vec3d Vector_normalize(vec3d &v) {
		float l = Vector_length(v);
		return { v.x / l , v.y / l, v.z / l, v.w };
	}


	bool OnUserCreate() override
	{

		meshCube.LoadFromObjectFile("mountains.obj");
		lookCamera = { 0.0f, 0.0f, 1.0f, 0.0f };
		lookCameraRight = { 1.0f, 0.0f, 0.0f, 0.0f };

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{



		mat4x4 pos = Matrix_MakeIdentity();
		pos.m[3][0] = 0.0f;
		pos.m[3][1] = 0.0f;
		pos.m[3][2] = 10.0f;

	
		if (GetKey(VK_UP).bHeld) {
			vCamera.y += 8.0f * fElapsedTime;
		}

		if (GetKey(VK_DOWN).bHeld) {
			vCamera.y-= 8.0f * fElapsedTime;
		}

		if (GetKey(VK_RIGHT).bHeld) {
			vCamera = Vector_sub(vCamera, Vector_mul(lookCameraRightRotated, 8.0f * fElapsedTime));
		}

		if (GetKey(VK_LEFT).bHeld) {
			vCamera = Vector_add(vCamera, Vector_mul(lookCameraRightRotated, 8.0f * fElapsedTime));
		}

		if (GetKey(L'A').bHeld) {
			fYaw -= 2.0f * fElapsedTime;
		}

		if (GetKey(L'D').bHeld) {
			fYaw += 2.0f * fElapsedTime;
		}


		mat4x4 matRotY = Matrix_MakeRotationY(fYaw);
		lookCameraRotated = MultiplyMatrixVector(lookCamera, matRotY);
		lookCameraRightRotated = MultiplyMatrixVector(lookCameraRight, matRotY);

		if (GetKey(L'S').bHeld)
			vCamera = Vector_sub(vCamera,Vector_mul(lookCameraRotated, 8.0f * fElapsedTime));

		if (GetKey(L'W').bHeld)
			vCamera = Vector_add(vCamera, Vector_mul(lookCameraRotated, 8.0f * fElapsedTime));

		// Clear Screen
		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);

		mesh TriToRaster;

		for (auto tri : meshCube.tris) {

			for (int i = 0; i < 3; i++) {
				//rotation 
				//position
				tri.p[i] = MultiplyMatrixVector(tri.p[i], pos);

			}
			vec3d line1 = Vector_sub(tri.p[1], tri.p[0]);
			vec3d line2 = Vector_sub(tri.p[2], tri.p[0]);
			vec3d normal = Vector_crossProduct(line1, line2);

			normal = Vector_normalize(normal);


			if ((Vector_dot(normal, Vector_sub(tri.p[0], vCamera))) < 0.0f) {



				float luminance;
				vec3d light_direction = { 0.0f, 1.0f, -1.0f, 0.0f };
				light_direction = Vector_normalize(light_direction);

				luminance = max(0.1f, (float)Vector_dot(light_direction, normal));


				tri.infos = GetColour(luminance);

				for (int i = 0; i < 3; i++) {
					tri.p[i] = Vector_toLocal(tri.p[i], vCamera, lookCameraRotated);
				}
				int nClippedTriangles = 0;
				triangle clipped[2];
				nClippedTriangles = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.2f }, { 0.0f, 0.0f, 1.0f }, tri, clipped[0], clipped[1]);

				for (int i = 0; i < nClippedTriangles; i++) {
					TriToRaster.tris.push_back(clipped[i]);
				}


			}
		}


		std::sort(TriToRaster.tris.begin(), TriToRaster.tris.end(), [](triangle v1, triangle v2) {
			float z1 = (v1.p[0].z + v1.p[1].z + v1.p[2].z) / 3.0f;
			float z2 = (v2.p[0].z + v2.p[1].z + v2.p[2].z) / 3.0f;

			return z2 < z1;

			});

		for (auto tri : TriToRaster.tris) {
			for (int i = 0; i < 3; i++) {


				if (tri.p[i].z != 0.0f) {
					tri.p[i] = Vector_div(tri.p[i], tri.p[i].z);
				}

				tri.p[i].x *= -1;
				tri.p[i].y *= -1;

				tri.p[i].x -= -1.0f;
				tri.p[i].y -= -1.0f;


				tri.p[i].y *= 0.5f * ScreenHeight();
				tri.p[i].x *=0.5f * ScreenWidth();
			}

			triangle clipped[2];
			list<triangle>listTriangles;
			listTriangles.push_back(tri);
			int nNewTriangles = 1;
			for (int p = 0; p < 4; p++) {

				int nTrisToAdd = 0;
				while (nNewTriangles > 0) {

					triangle test = listTriangles.front();
					listTriangles.pop_front();
					nNewTriangles--;
	

					switch (p)
					{
					case 0: nTrisToAdd = Triangle_ClipAgainstPlane({ 0, 0, 0, 0 }, { 0, 1.0f, 0, 0 }, test, clipped[0], clipped[1]); break;
					case 1: nTrisToAdd = Triangle_ClipAgainstPlane({ 0, (float)ScreenHeight() -1, 0, 0 }, { 0, -1.0f, 0, 0 }, test, clipped[0], clipped[1]); break;
					case 2: nTrisToAdd = Triangle_ClipAgainstPlane({ 0, 0, 0, 0 }, { 1.0f, 0, 0, 0 }, test, clipped[0], clipped[1]); break;
					case 3: nTrisToAdd = Triangle_ClipAgainstPlane({ (float)ScreenWidth() -1, 0, 0, 0 }, { -1.0f, 0, 0, 0 }, test, clipped[0], clipped[1]); break;
					}

					for (int s = 0; s < nTrisToAdd; s++) {
						listTriangles.push_back(clipped[s]);
					}
				}

				nNewTriangles = listTriangles.size();
			}
			for (auto triClipped : listTriangles) {

				FillTriangle(triClipped.p[0].x, triClipped.p[0].y,
					triClipped.p[1].x, triClipped.p[1].y,
					triClipped.p[2].x, triClipped.p[2].y,
					triClipped.infos.Char.UnicodeChar, triClipped.infos.Attributes);

				/*DrawTriangle(triClipped.p[0].x, triClipped.p[0].y,
					triClipped.p[1].x, triClipped.p[1].y,
					triClipped.p[2].x, triClipped.p[2].y,
					PIXEL_SOLID, FG_WHITE);*/
			}


		}


		return true;
	}


};




int main()
{
	olcEngine3D demo;
	if (demo.ConstructConsole(240, 240, 4, 4))
		demo.Start();

	return 0;
}
