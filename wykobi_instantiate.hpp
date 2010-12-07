/*
(***********************************************************************)
(*                                                                     *)
(* Wykobi Computational Geometry Library                               *)
(* Release Version 0.0.4                                               *)
(* http://www.wykobi.com                                               *)
(* Copyright (c) 2005-2009 Arash Partow, All Rights Reserved.          *)
(*                                                                     *)
(* The Wykobi computational geometry library and its components are    *)
(* supplied under the terms of the General Wykobi License agreement.   *)
(* The contents of the Wykobi computational geometry library and its   *)
(* components may not be copied or disclosed except in accordance with *)
(* the terms of that agreement.                                        *)
(*                                                                     *)
(* URL: http://www.wykobi.com/license.html                             *)
(*                                                                     *)
(***********************************************************************)
*/


#ifndef INCLUDE_WYKOBI_INSTANTIATE
#define INCLUDE_WYKOBI_INSTANTIATE

#include "wykobi.hpp"
#include "wykobi_nd.hpp"
#include "wykobi_utilities.hpp"

namespace wykobi
{

   #define INSTANTIATE_WYKOBI(T,InputIterator2d,InputIterator3d,OutputIterator2d,OutputIterator3d)\
      template inline int orientation<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& px, const T& py);\
      template inline int orientation<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& px, const T& py, const T& pz);\
      template inline int robust_orientation<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& px, const T& py);\
      template inline int robust_orientation<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& px, const T& py, const T& pz);\
      template inline int orientation<T>(const point2d<T>& point1, const point2d<T>& point2, const T& px, const T& py);\
      template inline int orientation<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline int orientation<T>(const line<T,2>& line, const point2d<T>& point);\
      template inline int orientation<T>(const segment<T,2>& segment, const point2d<T>& point);\
      template inline int orientation<T>(const triangle<T,2>& triangle);\
      template inline int orientation<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const T& px, const T& py, const T& pz);\
      template inline int orientation<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const point3d<T>& point4);\
      template inline int orientation<T>(const triangle<T,3>& triangle, const point3d<T>& point);\
      template inline bool differing_orientation<T>(const T& x1,  const T& y1, const T& x2,  const T& y2, const T& p1x, const T& p1y, const T& p2x, const T& p2y);\
      template inline bool differing_orientation<T>(const point2d<T>& p1, const point2d<T>& p2, const point2d<T>& q1, const point2d<T>& q2);\
      template inline int in_circle<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& px, const T& py);\
      template inline int in_circle<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4);\
      template inline int in_circle<T>(const triangle<T,2>& triangle, const point2d<T>& point);\
      template inline int in_sphere<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4, const T& px, const T& py, const T& pz);\
      template inline int in_sphere<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const point3d<T>& point4, const point3d<T>& point5);\
      template inline int in_sphere<T>(const quadix<T,3>& quadix, const point3d<T>& point);\
      template inline T signed_area<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& px, const T& py);\
      template inline T signed_area<T>(const point2d<T>& point1, const point2d<T>& point2, const T& px, const T& py);\
      template inline T signed_area<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline T signed_area<T>(const segment<T,2>& segment, const point2d<T>& point);\
      template inline T signed_volume<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& px, const T& py, const T& pz);\
      template inline T signed_volume<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const T& px, const T& py, const T& pz);\
      template inline T signed_volume<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const point3d<T>& point4);\
      template inline T signed_volume<T>(const triangle<T,3>& triangle, const point3d<T>& point);\
      template inline bool collinear<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& epsilon);\
      template inline bool collinear<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& epsilon);\
      template inline bool collinear<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline bool collinear<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3);\
      template inline bool robust_collinear<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& epsilon);\
      template inline bool robust_collinear<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const T& epsilon);\
      template inline bool robust_collinear<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& epsilon);\
      template inline bool robust_collinear<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const T& epsilon);\
      template inline bool robust_collinear<T>(const line<T,2>& line, const point2d<T>& point, const T& epsilon);\
      template inline bool robust_collinear<T>(const line<T,3>& line, const point3d<T>& point, const T& epsilon);\
      template inline bool is_point_collinear<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& px, const T& py, const bool robust);\
      template inline bool is_point_collinear<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const bool robust);\
      template inline bool is_point_collinear<T>(const point2d<T>& point1, const point2d<T>& point2, const T& px, const T& py, const bool robust);\
      template inline bool is_point_collinear<T>(const segment<T,2>& segment, const point2d<T>&point, const bool robust);\
      template inline bool is_point_collinear<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& px, const T& py, const T& pz, const bool robust);\
      template inline bool is_point_collinear<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const bool robust);\
      template inline bool is_point_collinear<T>(const segment<T,3>& segment, const point3d<T>&point, const bool robust);\
      template inline bool robust_coplanar<T>(const point3d<T> point1, const point3d<T> point2, const point3d<T> point3, const point3d<T> point4, const T& epsilon);\
      template inline bool coplanar<T>(const ray<T,3>& ray1, const ray<T,3>& ray2);\
      template inline bool coplanar<T>(const segment<T,3>& segment1, const segment<T,3>& segment2);\
      template inline bool coplanar<T>(const line<T,3>& line1, const line<T,3>& line2);\
      template inline bool coplanar<T>(const triangle<T,3>& triangle1, const triangle<T,3>& triangle2);\
      template inline bool coplanar<T>(const quadix<T,3>& quadix1, const quadix<T,3>& quadix2);\
      template inline bool cocircular<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, const T& epsilon);\
      template inline bool cocircular<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4, const T& epsilon);\
      template inline bool cocircular<T>(const triangle<T,2>& triangle, const point2d<T>& point, const T& epsilon);\
      template inline bool cocircular<T>(const circle<T>& circle, const point2d<T>& point, const T& epsilon);\
      template inline bool is_skinny_triangle<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline bool is_skinny_triangle<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline bool is_skinny_triangle<T>(const triangle<T,2>& triangle);\
      template inline bool intersect<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline bool intersect<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, T& ix,T& iy);\
      template inline bool intersect<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4);\
      template inline bool intersect<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4, point2d<T>& int_point);\
      template inline bool intersect<T>(const segment<T,2>& segment1, const segment<T,2>& segment2);\
      template inline bool intersect<T>(const segment<T,2>& segment1, const segment<T,2>& segment2,T& ix, T& iy);\
      template inline bool intersect<T>(const segment<T,2>& segment1, const segment<T,2>& segment2,point2d<T>& i_point);\
      template inline bool intersect<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4, const T& fuzzy);\
      template inline bool intersect<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const point3d<T>& point4, const T& fuzzy);\
      template inline bool intersect<T>(const segment<T,3>& segment1, const segment<T,3>& segment2, const T& fuzzy);\
      template inline bool intersect<T>(const segment<T,2>& segment, const rectangle<T>& rectangle);\
      template inline bool intersect<T>(const segment<T,2>& segment, const triangle<T,2>& triangle);\
      template inline bool intersect<T>(const segment<T,2>& segment, const quadix<T,2>& quadix);\
      template inline bool intersect<T>(const segment<T,2>& segment, const line<T,2>& line);\
      template inline bool intersect<T>(const segment<T,2>& segment, const circle<T>& circle);\
      template inline bool intersect<T>(const segment<T,2>& segment, const quadratic_bezier<T,2>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const segment<T,2>& segment, const cubic_bezier<T,2>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const segment<T,3>& segment, const line<T,3>& line, const T& fuzzy);\
      template inline bool intersect<T>(const segment<T,3>& segment, const box<T,3>& box);\
      template inline bool intersect<T>(const segment<T,3>& segment, const sphere<T>& sphere);\
      template inline bool intersect<T>(const segment<T,3>& segment, const plane<T,3>& plane);\
      template inline bool intersect<T>(const segment<T,3>& segment, const quadratic_bezier<T,3>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const segment<T,3>& segment, const cubic_bezier<T,3>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const line<T,2>& line, const triangle<T,2>& triangle);\
      template inline bool intersect<T>(const line<T,2>& line, const quadix<T,2>& quadix);\
      template inline bool intersect<T>(const line<T,2>& line1, const line<T,2>& line2);\
      template inline bool intersect<T>(const line<T,2>& line, const circle<T>& circle);\
      template inline bool intersect<T>(const line<T,2>& line, const quadratic_bezier<T,2>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const line<T,2>& line, const cubic_bezier<T,2>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const line<T,3>& line, const triangle<T,3>& triangle);\
      template inline bool intersect<T>(const line<T,3>& line, const plane<T,3>& plane);\
      template inline bool intersect<T>(const line<T,3>& line, const sphere<T>& sphere);\
      template inline bool intersect<T>(const line<T,3>& line, const quadratic_bezier<T,3>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const line<T,3>& line, const cubic_bezier<T,3>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const triangle<T,2>& triangle1, const triangle<T,2>& triangle2);\
      template inline bool intersect<T>(const triangle<T,2>& triangle, const circle<T>& circle);\
      template inline bool intersect<T>(const triangle<T,2>& triangle, const rectangle<T>& rectangle);\
      template inline bool intersect<T>(const triangle<T,2>& triangle, const quadratic_bezier<T,2>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const triangle<T,2>& triangle, const cubic_bezier<T,2>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const rectangle<T>& rectangle1, const rectangle<T>& rectangle2);\
      template inline bool intersect<T>(const rectangle<T>& rectangle, const quadratic_bezier<T,2>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const rectangle<T>& rectangle, const cubic_bezier<T,2>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const rectangle<T>& rectangle, const circle<T>& circle);\
      template inline bool intersect<T>(const quadix<T,2>& quadix, const quadratic_bezier<T,2>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const quadix<T,2>& quadix, const cubic_bezier<T,2>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const circle<T>& circle1, const circle<T>& circle2);\
      template inline bool intersect<T>(const circle<T>& circle, const quadratic_bezier<T,2>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const circle<T>& circle, const cubic_bezier<T,2>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const box<T,3>& box, const sphere<T>& sphere);\
      template inline bool intersect<T>(const sphere<T>& sphere1, const sphere<T>& sphere2);\
      template inline bool intersect<T>(const sphere<T>& sphere, const quadratic_bezier<T,3>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const sphere<T>& sphere, const cubic_bezier<T,3>& bezier, const std::size_t& steps);\
      template inline bool intersect<T>(const ray<T,2>& ray1, const ray<T,2>& ray2);\
      template inline bool intersect<T>(const ray<T,3>& ray1, const ray<T,3>& ray2);\
      template inline bool intersect<T>(const ray<T,2>& ray, const segment<T,2>& segment);\
      template inline bool intersect<T>(const ray<T,3>& ray, const segment<T,3>& segment);\
      template inline bool intersect<T>(const ray<T,2>& ray, const rectangle<T>& rectangle);\
      template inline bool intersect<T>(const ray<T,3>& ray, const box<T,3>& box);\
      template inline bool intersect<T>(const ray<T,2>& ray, const triangle<T,2>& triangle);\
      template inline bool intersect<T>(const ray<T,3>& ray, const triangle<T,3>& triangle);\
      template inline bool intersect<T>(const ray<T,2>& ray, const quadix<T,2>& quadix);\
      template inline bool intersect<T>(const ray<T,2>& ray, const circle<T>& circle);\
      template inline bool intersect<T>(const ray<T,3>& ray, const sphere<T>& sphere);\
      template inline bool intersect<T>(const ray<T,3>& ray, const plane<T,3>& plane);\
      template inline bool intersect<T>(const ray<T,2>& ray, const polygon<T,2>& polygon);\
      template inline bool intersect<T>(const plane<T,3>& plane1, const plane<T,3>& plane2);\
      template inline bool intersect<T>(const plane<T,3>& plane, const sphere<T>& sphere);\
      template inline bool intersect<T>(const plane<T,3>& plane, const line<T,3>& line);\
      template inline bool simple_intersect<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline bool simple_intersect<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4);\
      template inline bool simple_intersect<T>(const segment<T,2>& segment1, const segment<T,2>& segment2);\
      template inline bool intersect_vertical_horizontal<T>(const segment<T,2>& segment1, const segment<T,2>& segment2);\
      template inline bool intersect_vertical_vertical<T>(const segment<T,2>& segment1, const segment<T,2>& segment2);\
      template inline bool intersect_horizontal_horizontal<T>(const segment<T,2>& segment1, const segment<T,2>& segment2);\
      template inline void intersection_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, T& ix, T& iy);\
      template inline void intersection_point<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4, T& ix, T& iy);\
      template inline point2d<T> intersection_point<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4);\
      template inline point2d<T> intersection_point<T>(const segment<T,2>& segment1, const segment<T,2>& segment2);\
      template inline void intersection_point<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4, T& ix, T& iy, T& iz, const T& fuzzy);\
      template inline void intersection_point<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const point3d<T>& point4, T& ix, T& iy, T& iz, const T& fuzzy);\
      template inline point3d<T> intersection_point<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const point3d<T>& point4, const T& fuzzy);\
      template inline point3d<T> intersection_point<T>(const segment<T,3>& segment1, const segment<T,3>& segment2, const T& fuzzy);\
      template inline point2d<T> intersection_point<T>(const segment<T,2>& segment, const line<T,2>& line);\
      template inline point3d<T> intersection_point<T>(const segment<T,3>& segment, const line<T,3>& line, const T& fuzzy);\
      template inline point3d<T> intersection_point<T>(const segment<T,3>& segment, const plane<T,3>& plane);\
      template inline void intersection_point<T,OutputIterator3d>(const segment<T,3>& segment, const sphere<T>& sphere, OutputIterator3d out);\
      template inline void intersection_point<T,OutputIterator2d>(const segment<T,2>& segment, const quadratic_bezier<T,2>& bezier, OutputIterator2d out, const std::size_t& steps);\
      template inline void intersection_point<T,OutputIterator2d>(const segment<T,2>& segment, const cubic_bezier<T,2>& bezier, OutputIterator2d out, const std::size_t& steps);\
      template inline point2d<T> intersection_point<T>(const line<T,2>& line1, const line<T,2>& line2);\
      template inline point3d<T> intersection_point<T>(const line<T,3>& line1, const line<T,3>& line2, const T& fuzzy);\
      template inline void intersection_point<T>(const circle<T>& circle1, const circle<T>& circle2, point2d<T>& point1, point2d<T>& point2);\
      template inline void intersection_point<T,OutputIterator2d>(const segment<T,2>& segment, const triangle<T,2>& triangle, OutputIterator2d out);\
      template inline void intersection_point<T>(const line<T,3>& line, const triangle<T,3>& triangle, point3d<T>& ipoint);\
      template inline point3d<T> intersection_point<T>(const line<T,3>& line, const plane<T,3>& plane);\
      template inline void intersection_point<T,OutputIterator2d>(const T& x1, const T& y1, const T& x2, const T& y2, const T& cx, const T& cy, const T& radius, OutputIterator2d out);\
      template inline void intersection_point<T,OutputIterator2d>(const segment<T,2>& segment, const circle<T>& circle, OutputIterator2d out);\
      template inline void intersection_point<T,OutputIterator2d>(const line<T,2>& line, const circle<T>& circle, OutputIterator2d out);\
      template inline void intersection_point<T,OutputIterator3d>(const line<T,3>& line, const sphere<T>& sphere, OutputIterator3d out);\
      template inline point2d<T> intersection_point<T>(const ray<T,2>& ray1, const ray<T,2>& ray2);\
      template inline point3d<T> intersection_point<T>(const ray<T,3>& ray, const triangle<T,3>& triangle);\
      template inline point3d<T> intersection_point<T>(const ray<T,3>& ray, const plane<T,3>& plane);\
      template inline void intersection_point<T,OutputIterator2d>(const ray<T,2>& ray, const circle<T>& circle, OutputIterator2d out);\
      template inline void intersection_point<T,OutputIterator3d>(const ray<T,3>& ray, const sphere<T>& sphere, OutputIterator3d out);\
      template inline void intersection_point_line_to_line<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4, T& ix, T& iy, T& iz, const T& fuzzy);\
      template inline T normalize_angle<T>(const T& angle);\
      template inline T vertical_mirror<T>(const T& angle);\
      template inline T horizontal_mirror<T>(const T& angle);\
      template inline unsigned int quadrant<T>(const T& angle);\
      template inline unsigned int quadrant<T>(const T& x, const T& y);\
      template inline unsigned int quadrant<T>(const point2d<T>& point);\
      template inline T vertex_angle<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline T vertex_angle<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline T vertex_angle<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3);\
      template inline T vertex_angle<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3);\
      template inline T oriented_vertex_angle<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const int orient);\
      template inline T oriented_vertex_angle<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const int orient);\
      template inline T cartesian_angle<T>(const T& x, const T& y);\
      template inline T cartesian_angle<T>(const point2d<T>& point);\
      template inline T robust_cartesian_angle<T>(const T& x, const T& y);\
      template inline T robust_cartesian_angle<T>(const point2d<T>& point);\
      template inline T cartesian_angle<T>(const T& x, const T& y, const T& ox, const T& oy);\
      template inline T cartesian_angle<T>(const point2d<T>& point, const point2d<T>& origin);\
      template inline T robust_cartesian_angle<T>(const T& x, const T& y, const T& ox, const T& oy);\
      template inline T robust_cartesian_angle<T>(const point2d<T>& point, const point2d<T>& origin);\
      template inline bool parallel<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, const T& epsilon);\
      template inline bool parallel<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4, const T& epsilon);\
      template inline bool parallel<T>(const segment<T,2>& segment1, const segment<T,2>& segment2, const T& epsilon);\
      template inline bool parallel<T>(const line<T,2>& line1, const line<T,2>& line2, const T& epsilon);\
      template inline bool parallel<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4, const T& epsilon);\
      template inline bool parallel<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const point3d<T>& point4, const T& epsilon);\
      template inline bool parallel<T>(const segment<T,3>& segment1, const segment<T,3>& segment2, const T& epsilon);\
      template inline bool parallel<T>(const line<T,3>& line1, const line<T,3>& line2, const T& epsilon);\
      template inline bool robust_parallel<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, const T& epsilon);\
      template inline bool robust_parallel<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4, const T& epsilon);\
      template inline bool robust_parallel<T>(const segment<T,2>& segment1, const segment<T,2>& segment2, const T& epsilon);\
      template inline bool robust_parallel<T>(const line<T,2>& line1, const line<T,2>& line2, const T& epsilon);\
      template inline bool robust_parallel<T>(const line<T,2>& line, const segment<T,2>& segment, const T& epsilon);\
      template inline bool robust_parallel<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4, const T& epsilon);\
      template inline bool robust_parallel<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const point3d<T>& point4, const T& epsilon);\
      template inline bool robust_parallel<T>(const segment<T,3>& segment1, const segment<T,3>& segment2, const T& epsilon);\
      template inline bool robust_parallel<T>(const line<T,3>& line1, const line<T,3>& line2, const T& epsilon);\
      template inline bool robust_parallel<T>(const line<T,3>& line, const segment<T,3>& segment, const T& epsilon);\
      template inline bool perpendicular<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, const T& epsilon);\
      template inline bool perpendicular<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4, const T& epsilon);\
      template inline bool perpendicular<T>(const segment<T,2>& segment1, const segment<T,2>& segment2, const T& epsilon);\
      template inline bool perpendicular<T>(const line<T,2>& line1, const line<T,2>& line2, const T& epsilon);\
      template inline bool perpendicular<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4, const T& epsilon);\
      template inline bool perpendicular<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const point3d<T>& point4, const T& epsilon);\
      template inline bool perpendicular<T>(const segment<T,3>& segment1, const segment<T,3>& segment2, const T& epsilon);\
      template inline bool perpendicular<T>(const line<T,3>& line1, const line<T,3>& line2, const T& epsilon);\
      template inline bool perpendicular<T>(const line<T,2>& line, const segment<T,2>& segment, const T& epsilon);\
      template inline bool robust_perpendicular<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, const T& epsilon);\
      template inline bool robust_perpendicular<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4, const T& epsilon);\
      template inline bool robust_perpendicular<T>(const segment<T,2>& segment1, const segment<T,2>& segment2, const T& epsilon);\
      template inline bool robust_perpendicular<T>(const line<T,2>& line1, const line<T,2>& line2, const T& epsilon);\
      template inline bool robust_perpendicular<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4, const T& epsilon);\
      template inline bool robust_perpendicular<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const point3d<T>& point4, const T& epsilon);\
      template inline bool robust_perpendicular<T>(const segment<T,3>& segment1, const segment<T,3>& segment2, const T& epsilon);\
      template inline bool robust_perpendicular<T>(const line<T,3>& line1, const line<T,3>& line2, const T& epsilon);\
      template inline bool robust_perpendicular<T>(const line<T,2>& line, const segment<T,2>& segment, const T& epsilon);\
      template inline bool line_to_line_intersect<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline bool line_to_line_intersect<T>(const line<T,2>& line1, const line<T,2>& line2);\
      template inline bool rectangle_to_rectangle_intersect<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline bool rectangle_to_rectangle_intersect<T>(const rectangle<T>& rectangle1, const rectangle<T>& rectangle2);\
      template inline bool box_to_box_intersect<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4);\
      template inline bool box_to_box_intersect<T>(const box<T,3>& box1, const box<T,3>& box2);\
      template inline bool rectangle_within_rectangle<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline bool rectangle_within_rectangle<T>(const rectangle<T>& rectangle1, const rectangle<T>& rectangle2);\
      template inline bool box_within_box<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4);\
      template inline bool box_within_box<T>(const box<T,3>& box1, const box<T,3>& box2);\
      template inline bool circle_within_rectangle<T>(const T& x, const T& y, const T& radius, const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline bool circle_within_rectangle<T>(const circle<T>& circle, const rectangle<T>& rectangle);\
      template inline bool triangle_within_rectangle<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, const T& x5, const T& y5);\
      template inline bool triangle_within_rectangle<T>(const triangle<T,2>& triangle, const rectangle<T>& rectangle);\
      template inline bool segment_within_rectangle<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline bool segment_within_rectangle<T>(const segment<T,2>& segment, const rectangle<T>& rectangle);\
      template inline bool quadix_within_rectangle<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, const T& x5, const T& y5, const T& x6, const T& y6);\
      template inline bool quadix_within_rectangle<T>(const quadix<T,2>& quadix, const rectangle<T>& rectangle);\
      template inline bool polygon_within_rectangle<T>(const polygon<T,2>& polygon, const rectangle<T>& rectangle);\
      template inline bool sphere_within_box<T>(const T& x, const T& y, const T& z, const T& radius, const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline bool sphere_within_box<T>(const sphere<T>& sphere, const box<T,3>& box);\
      template inline bool triangle_within_box<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4, const T& x5, const T& y5, const T& z5);\
      template inline bool triangle_within_box<T>(const triangle<T,3>& triangle, const box<T,3>& box);\
      template inline bool segment_within_box<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4) ;\
      template inline bool segment_within_box<T>(const segment<T,3>& segment, const box<T,3>& box);\
      template inline bool quadix_within_box<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4, const T& x5, const T& y5, const T& z5, const T& x6, const T& y6, const T& z6);\
      template inline bool quadix_within_box<T>(const quadix<T,3>& quadix, const box<T,3>& box);\
      template inline bool polygon_within_box<T>(const polygon<T,3>& polygon, const box<T,3>& box);\
      template inline bool circle_in_circle<T>(const circle<T>& circle1, const circle<T>& circle2);\
      template inline bool is_tangent<T>(const segment<T,2>& segment, const circle<T>& circle);\
      template inline bool point_of_reflection<T>(const T& sx1, const T& sy1, const T& sx2, const T& sy2, const T& p1x, const T& p1y, const T& p2x, const T& p2y, T& rpx, T& rpy);\
      template inline bool point_of_reflection<T>(const segment<T,2>& segment, const point2d<T>&point1, const point2d<T>&point2, point2d<T>&reflection_point);\
      template inline segment<T,2> edge<T>(const triangle<T,2>& triangle, const std::size_t& edge_index);\
      template inline segment<T,3> edge<T>(const triangle<T,3>& triangle, const std::size_t& edge_index);\
      template inline segment<T,2> edge<T>(const quadix<T,2>& quadix, const std::size_t& edge_index);\
      template inline segment<T,3> edge<T>(const quadix<T,3>& quadix, const std::size_t& edge_index);\
      template inline segment<T,2> edge<T>(const rectangle<T>& rectangle, const std::size_t& edge);\
      template inline segment<T,2> edge<T>(const polygon<T,2>& polygon, const std::size_t& edge);\
      template inline segment<T,3> edge<T>(const polygon<T,3>& polygon, const std::size_t& edge);\
      template inline segment<T,2> opposing_edge<T>(const triangle<T,2>& triangle, const std::size_t& corner);\
      template inline segment<T,3> opposing_edge<T>(const triangle<T,3>& triangle, const std::size_t& corner);\
      template inline segment<T,2> reverse_segment<T>(const segment<T,2>& segment);\
      template inline segment<T,3> reverse_segment<T>(const segment<T,3>& segment);\
      template inline point2d<T> rectangle_corner<T>(const rectangle<T>& rectangle, const std::size_t& corner_index);\
      template inline point3d<T> box_corner<T>(const box<T,3>& box, const std::size_t& corner_index);\
      template inline line<T,2> triangle_bisector<T>(const triangle<T,2>& triangle, const std::size_t& bisector);\
      template inline line<T,3> triangle_bisector<T>(const triangle<T,3>& triangle, const std::size_t& bisector);\
      template inline line<T,2> triangle_external_bisector<T>(const triangle<T,2>& triangle, const std::size_t& corner, const std::size_t& opposing_corner);\
      template inline line<T,3> triangle_external_bisector<T>(const triangle<T,3>& triangle, const std::size_t& corner, const std::size_t& opposing_corner);\
      template inline line<T,2> triangle_median<T>(const triangle<T,2>& triangle, const std::size_t& median);\
      template inline line<T,3> triangle_median<T>(const triangle<T,3>& triangle, const std::size_t& median);\
      template inline line<T,2> triangle_symmedian<T>(const triangle<T,2>& triangle, const std::size_t& symmedian);\
      template inline line<T,3> triangle_symmedian<T>(const triangle<T,3>& triangle, const std::size_t& symmedian);\
      template inline line<T,2> euler_line<T>(const triangle<T,2>& triangle);\
      template inline line<T,3> euler_line<T>(const triangle<T,3>& triangle);\
      template inline point2d<T> exmedian_point<T>(const triangle<T,2>& triangle, const std::size_t& corner);\
      template inline point3d<T> exmedian_point<T>(const triangle<T,3>& triangle, const std::size_t& corner);\
      template inline point2d<T> feuerbach_point<T>(const triangle<T,2>& triangle);\
      template inline line<T,2> confined_triangle_median<T>(const triangle<T,2>& triangle,const point2d<T>& point, const std::size_t& median);\
      template inline line<T,3> confined_triangle_median<T>(const triangle<T,3>& triangle,const point3d<T>& point, const std::size_t& median);\
      template inline line<T,2> create_parallel_line_on_point<T>(const line<T,2>& line, const point2d<T>& point);\
      template inline line<T,3> create_parallel_line_on_point<T>(const line<T,3>& line, const point3d<T>& point);\
      template inline segment<T,2> create_parallel_segment_on_point<T>(const line<T,2>& line, const point2d<T>& point);\
      template inline segment<T,3> create_parallel_segment_on_point<T>(const line<T,3>& line, const point3d<T>& point);\
      template inline bool point_in_rectangle<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline bool point_in_rectangle<T>(const point2d<T>& point, const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline bool point_in_rectangle<T>(const T& px, const T& py, const rectangle<T>& rectangle);\
      template inline bool point_in_rectangle<T>(const point2d<T>& point, const rectangle<T>& rectangle);\
      template inline bool point_in_rectangle<T>(const point2d<T>& point, const point2d<T>& rect_point1, point2d<T>& rect_point2);\
      template inline bool point_in_rectangle<T>(const point2d<T>& point, const segment<T,2>& segment);\
      template inline bool point_in_box<T>(const T& px, const T& py, const T& pz, const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline bool point_in_box<T>(const point3d<T>& point, const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline bool point_in_box<T>(const T& px, const T& py, const T& pz, const box<T,3>& box);\
      template inline bool point_in_box<T>(const point3d<T>& point, const box<T,3>& box);\
      template inline bool point_in_box<T>(const point3d<T>& point, const point3d<T>& box_point1, const point3d<T>& box_point2);\
      template inline bool point_in_box<T>(const point3d<T>& point, const segment<T,3>& segment);\
      template inline bool point_in_triangle<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline bool point_in_triangle<T>(const point2d<T>& point, const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline bool point_in_triangle<T>(const T& px, const T& py, const triangle<T,2>& triangle);\
      template inline bool point_in_triangle<T>(const point2d<T>& point, const triangle<T,2>& triangle);\
      template inline bool point_in_quadix<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline bool point_in_quadix<T>(const point2d<T>& point, const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4);\
      template inline bool point_in_quadix<T>(const T& px, const T& py, const quadix<T,2>& quadix);\
      template inline bool point_in_quadix<T>(const point2d<T>& point, const quadix<T,2>& quadix);\
      template inline bool point_in_circle<T>(const T& px, const T& py, const T& cx, const T& cy, const T& radius);\
      template inline bool point_in_circle<T>(const T& px, const T& py, const circle<T>& circle);\
      template inline bool point_in_circle<T>(const point2d<T>& point, const circle<T>& circle);\
      template inline bool point_in_sphere<T>(const T& px, const T& py, const T& pz, const T& cx, const T& cy, const T& cz, const T& radius);\
      template inline bool point_in_sphere<T>(const T& px, const T& py, const T& pz, const sphere<T>& sphere);\
      template inline bool point_in_sphere<T>(const point3d<T>& point, const sphere<T>& sphere);\
      template inline bool point_in_three_point_circle<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline bool point_in_three_point_circle<T>(const point2d<T>& point, const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline bool point_in_three_point_circle<T>(const point2d<T>& point, const triangle<T,2> triangle);\
      template inline bool point_in_focus_area<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline bool point_in_focus_area<T>(const point2d<T>& point, const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline bool point_on_segment<T>(const point2d<T>& point, const segment<T,2>& segment);\
      template inline bool point_on_segment<T>(const point3d<T>& point, const segment<T,3>& segment);\
      template inline bool point_on_ray<T>(const T& px, const T& py, const T& ox, const T& oy, const T& dx, const T& dy);\
      template inline bool point_on_ray<T>(const T& px, const T& py, const T& pz, const T& ox, const T& oy, const T& oz, const T& dx, const T& dy, const T& dz);\
      template inline bool point_on_ray<T>(const point2d<T>& point, const ray<T,2>& ray);\
      template inline bool point_on_ray<T>(const point3d<T>& point, const ray<T,3>& ray);\
      template inline bool point_on_rectangle<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline bool point_on_rectangle<T>(const point2d<T>& point, const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline bool point_on_rectangle<T>(const T& px, const T& py, const rectangle<T>& rectangle);\
      template inline bool point_on_rectangle<T>(const point2d<T>& point, const rectangle<T>& rectangle);\
      template inline bool point_on_triangle<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline bool point_on_triangle<T>(const point2d<T>& point, const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline bool point_on_triangle<T>(const T& px, const T& py, const triangle<T,2>& triangle);\
      template inline bool point_on_triangle<T>(const point2d<T>& point, const triangle<T,2>& triangle);\
      template inline bool point_on_quadix<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline bool point_on_quadix<T>(const point2d<T>& point, const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4);\
      template inline bool point_on_quadix<T>(const T& px, const T& py, const quadix<T,2>& quadix);\
      template inline bool point_on_quadix<T>(const point2d<T>& point, const quadix<T,2>& quadix);\
      template inline bool point_on_circle<T>(const T& px, const T& py, const T& cx, const T& cy, const T& radius);\
      template inline bool point_on_circle<T>(const T& px, const T& py, const circle<T>& circle);\
      template inline bool point_on_circle<T>(const point2d<T>& point, const circle<T>& circle);\
      template inline bool point_on_bezier<T>(const point2d<T>& point, const quadratic_bezier<T,2>& bezier, const std::size_t& steps, const T& fuzzy);\
      template inline bool point_on_bezier<T>(const point2d<T>& point, const cubic_bezier<T,2>& bezier, const std::size_t& steps, const T& fuzzy);\
      template inline bool point_on_bezier<T>(const point3d<T>& point, const quadratic_bezier<T,3>& bezier, const std::size_t& steps, const T& fuzzy);\
      template inline bool point_on_bezier<T>(const point3d<T>& point, const cubic_bezier<T,3>& bezier, const std::size_t& steps, const T& fuzzy);\
      template inline point2d<T> isogonal_conjugate<T>(const point2d<T>& point, const triangle<T,2>& triangle);\
      template inline point3d<T> isogonal_conjugate<T>(const point3d<T>& point, const triangle<T,3>& triangle);\
      template inline point2d<T> cyclocevian_conjugate<T>(const point2d<T>& point, const triangle<T,2>& triangle);\
      template inline point2d<T> symmedian_point<T>(const triangle<T,2>& triangle);\
      template inline point3d<T> symmedian_point<T>(const triangle<T,3>& triangle);\
      template inline void create_equilateral_triangle<T>(const T& x1, const T& y1, const T& x2, const T& y2, T& x3, T& y3);\
      template inline void create_equilateral_triangle<T>(const point2d<T>& point1, const point2d<T>& point2, point2d<T>& point3);\
      template inline triangle<T,2> create_equilateral_triangle<T>(const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline triangle<T,2> create_equilateral_triangle<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline triangle<T,2> create_equilateral_triangle<T>(const T& cx, const T& cy, const T& side_length);\
      template inline triangle<T,2> create_equilateral_triangle<T>(const point2d<T>& center_point, const T& side_length);\
      template inline triangle<T,2> create_isosceles_triangle<T>(const point2d<T>& point1, const point2d<T>& point2, const T& angle);\
      template inline triangle<T,2> create_isosceles_triangle<T>(const segment<T,2>& segment, const T& angle);\
      template inline triangle<T,2> create_triangle<T>(const point2d<T>& point1, const point2d<T>& point2, const T& angle1, const T& angle2);\
      template inline triangle<T,2> create_triangle<T>(const segment<T,2>& segment, const T& angle1, const T& angle2);\
      template inline triangle<T,2> create_morley_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,2> create_cevian_triangle<T>(const triangle<T,2>& triangle, const point2d<T>& point);\
      template inline triangle<T,3> create_cevian_triangle<T>(const triangle<T,3>& triangle, const point3d<T>& point);\
      template inline triangle<T,2> create_anticevian_triangle<T>(const triangle<T,2>& triangle, const point2d<T>& point);\
      template inline triangle<T,3> create_anticevian_triangle<T>(const triangle<T,3>& triangle, const point3d<T>& point);\
      template inline triangle<T,2> create_anticomplementary_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,3> create_anticomplementary_triangle<T>(const triangle<T,3>& triangle);\
      template inline triangle<T,2> create_inner_napoleon_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,2> create_outer_napoleon_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,2> create_inner_vecten_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,2> create_outer_vecten_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,2> create_medial_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,3> create_medial_triangle<T>(const triangle<T,3>& triangle);\
      template inline triangle<T,2> create_contact_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,3> create_contact_triangle<T>(const triangle<T,3>& triangle);\
      template inline triangle<T,2> create_symmedial_triangle<T>(const triangle<T,2>& triangle, const point2d<T>& point);\
      template inline triangle<T,2> create_orthic_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,3> create_orthic_triangle<T>(const triangle<T,3>& triangle);\
      template inline triangle<T,2> create_pedal_triangle<T>(const point2d<T>& point, const triangle<T,2>& triangle);\
      template inline triangle<T,3> create_pedal_triangle<T>(const point3d<T>& point, const triangle<T,3>& triangle);\
      template inline triangle<T,2> create_antipedal_triangle<T>(const point2d<T>& point, const triangle<T,2>& triangle);\
      template inline triangle<T,2> create_excentral_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,3> create_excentral_triangle<T>(const triangle<T,3>& triangle);\
      template inline triangle<T,2> create_incentral_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,3> create_incentral_triangle<T>(const triangle<T,3>& triangle);\
      template inline triangle<T,2> create_extouch_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,3> create_extouch_triangle<T>(const triangle<T,3>& triangle);\
      template inline triangle<T,2> create_feuerbach_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,2> create_circumcevian_triangle<T>(const triangle<T,2>& triangle, const point2d<T>& point);\
      template inline triangle<T,2> create_circummedial_triangle<T>(const triangle<T,2>& triangle);\
      template inline triangle<T,2> create_first_brocard_triangle<T>(const triangle<T,2>& triangle);\
      template inline void create_right_triangle<T>(const wykobi::point2d<T>& p1, const wykobi::point2d<T>& p2, wykobi::point2d<T>& c1, wykobi::point2d<T>& c2);\
      template inline void create_equilateral_quadix<T>(const T& x1, const T& y1, const T& x2, const T& y2, T& x3, T& y3, T& x4, T& y4);\
      template inline void create_equilateral_quadix<T>(const point2d<T>& point1, const point2d<T>& point2, point2d<T>& point3, point2d<T>& point4);\
      template inline quadix<T,2> create_equilateral_quadix<T>(const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline quadix<T,2> create_equilateral_quadix<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline quadix<T,2> create_equilateral_quadix<T>(const segment<T,2>& segment);\
      template inline quadix<T,2> create_equilateral_quadix<T>(const T& cx, const T& cy, const T& side_length);\
      template inline quadix<T,2> create_equilateral_quadix<T>(const point2d<T>& center_point, const T& side_length);\
      template inline void torricelli_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, T& px, T& py);\
      template inline point2d<T> torricelli_point<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline point2d<T> torricelli_point<T>(const triangle<T,2>& triangle);\
      template inline void incenter<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, T& px, T& py);\
      template inline void incenter<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, T& px, T& py, T& pz);\
      template inline point2d<T> incenter<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline point3d<T> incenter<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3);\
      template inline point2d<T> incenter<T>(const triangle<T,2>& triangle);\
      template inline point3d<T> incenter<T>(const triangle<T,3>& triangle);\
      template inline void circumcenter<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, T& px, T& py);\
      template inline void circumcenter<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, T& px, T& py, T& pz);\
      template inline point2d<T> circumcenter<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline point3d<T> circumcenter<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3);\
      template inline point2d<T> circumcenter<T>(const triangle<T,2>& triangle);\
      template inline point3d<T> circumcenter<T>(const triangle<T,3>& triangle);\
      template inline circle<T> circumcircle<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline circle<T> circumcircle<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline circle<T> circumcircle<T>(const triangle<T,2>& triangle);\
      template inline sphere<T> circumsphere<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3);\
      template inline sphere<T> circumsphere<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3);\
      template inline sphere<T> circumsphere<T>(const triangle<T,3>& triangle);\
      template inline circle<T> inscribed_circle<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline circle<T> inscribed_circle<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline circle<T> inscribed_circle<T>(const triangle<T,2>& triangle);\
      template inline sphere<T> inscribed_sphere<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3);\
      template inline sphere<T> inscribed_sphere<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3);\
      template inline sphere<T> inscribed_sphere<T>(const triangle<T,3>& triangle);\
      template inline circle<T> nine_point_circle<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline circle<T> nine_point_circle<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline circle<T> nine_point_circle<T>(const triangle<T,2>& triangle);\
      template inline point2d<T> orthocenter<T>(const triangle<T,2>& triangle);\
      template inline point3d<T> orthocenter<T>(const triangle<T,3>& triangle);\
      template inline point2d<T> excenter<T>(const triangle<T,2>& triangle, const std::size_t& corner);\
      template inline point3d<T> excenter<T>(const triangle<T,3>& triangle, const std::size_t& corner);\
      template inline circle<T> excircle<T>(const triangle<T,2>& triangle, const std::size_t& corner);\
      template inline circle<T> mandart_circle<T>(const triangle<T,2>& triangle);\
      template inline circle<T> brocard_circle<T>(const triangle<T,2>& triangle);\
      template inline circle<T> invert_circle_across_circle<T>(const circle<T>& circle1, const circle<T>& circle2);\
      template inline sphere<T> invert_sphere_across_sphere<T>(const sphere<T>& sphere1, const sphere<T>& sphere2);\
      template inline void circle_tangent_points<T>(const circle<T>& circle, const point2d<T>& point, point2d<T>& point1, point2d<T>& point2);\
      template inline line<T,2> tangent_line<T>(const circle<T>& circle, const point2d<T>& point);\
      template inline line<T,2> create_line_from_bisector<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline segment<T,2> create_segment_from_bisector<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline line<T,3> create_line_from_bisector<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3);\
      template inline segment<T,3> create_segment_from_bisector<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3);\
      template inline line<T,2> create_line_from_bisector<T>(const point2d<T>& point1,const point2d<T>& point2,const point2d<T>& point3);\
      template inline segment<T,2> create_segment_from_bisector<T>(const point2d<T>& point1,const point2d<T>& point2,const point2d<T>& point3);\
      template inline ray<T,2> create_ray_from_bisector<T>(const point2d<T>& point1,const point2d<T>& point2,const point2d<T>& point3);\
      template inline line<T,3> create_line_from_bisector<T>(const point3d<T>& point1,const point3d<T>& point2,const point3d<T>& point3);\
      template inline segment<T,3> create_segment_from_bisector<T>(const point3d<T>& point1,const point3d<T>& point2,const point3d<T>& point3);\
      template inline ray<T,3> create_ray_from_bisector<T>(const point3d<T>& point1,const point3d<T>& point2,const point3d<T>& point3);\
      template inline line<T,2> create_perpendicular_bisector<T>(const T& x1, const T& y1,const T& x2, const T& y2);\
      template inline line<T,2> create_perpendicular_bisector<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline line<T,2> create_perpendicular_bisector<T>(const segment<T,2>& segment);\
      template inline line<T,2> create_perpendicular_line_at_end_point<T>(const line<T,2>& line);\
      template inline void closest_point_on_segment_from_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& px, const T& py, T& nx, T& ny);\
      template inline void closest_point_on_segment_from_point<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& px, const T& py, const T& pz, T& nx, T& ny, T& nz);\
      template inline void closest_point_on_line_from_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& px, const T& py, T& nx, T& ny);\
      template inline void closest_point_on_line_from_point<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& px, const T& py, const T& pz, T& nx, T& ny, T& nz);\
      template inline void order_sensitive_closest_point_on_segment_from_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& px, const T& py, T& nx, T& ny);\
      template inline void order_sensitive_closest_point_on_segment_from_point<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& px, const T& py, const T& pz, T& nx, T& ny, T& nz);\
      template inline void order_sensitive_closest_point_on_line_from_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& px, const T& py, T& nx, T& ny);\
      template inline void order_sensitive_closest_point_on_line_from_point<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& px, const T& py, const T& pz, T& nx, T& ny, T& nz);\
      template inline void closest_point_on_ray_from_point<T>(const T& ox, const T& oy, const T& dx, const T& dy, const T& px, const T& py, T& nx, T& ny);\
      template inline void closest_point_on_ray_from_point<T>(const T& ox, const T& oy, const T& oz, const T& dx, const T& dy, const T& dz, const T& px, const T& py, const T& pz, T& nx, T& ny, T& nz);\
      template inline point2d<T> closest_point_on_segment_from_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& px, const T& py);\
      template inline point3d<T> closest_point_on_segment_from_point<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& px, const T& py, const T& pz);\
      template inline point2d<T> closest_point_on_segment_from_point<T>(const segment<T,2>& segment, const point2d<T>& point);\
      template inline point3d<T> closest_point_on_segment_from_point<T>(const segment<T,3>& segment, const point3d<T>& point);\
      template inline point2d<T> closest_point_on_line_from_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& px, const T& py);\
      template inline point3d<T> closest_point_on_line_from_point<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& px, const T& py, const T& pz);\
      template inline point2d<T> closest_point_on_line_from_point<T>(const line<T,2>& line, const point2d<T>& point);\
      template inline point3d<T> closest_point_on_line_from_point<T>(const line<T,3>& line, const point3d<T>& point);\
      template inline point2d<T> closest_point_on_ray_from_point<T>(const T& ox, const T& oy, const T& dx, const T& dy, const T& px, const T& py);\
      template inline point3d<T> closest_point_on_ray_from_point<T>(const T& ox, const T& oy, const T& oz, const T& dx, const T& dy, const T& dz, const T& px, const T& py, const T& pz);\
      template inline point2d<T> closest_point_on_ray_from_point<T>(const ray<T,2>& ray, const point2d<T>& point);\
      template inline point3d<T> closest_point_on_ray_from_point<T>(const ray<T,3>& ray, const point3d<T>& point);\
      template inline void closest_point_on_triangle_from_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& px, const T& py, T& nx, T& ny);\
      template inline point2d<T> closest_point_on_triangle_from_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& px, const T& py);\
      template inline point2d<T> closest_point_on_triangle_from_point<T>(const triangle<T,2>& triangle, const T& px, const T& py);\
      template inline point2d<T> closest_point_on_triangle_from_point<T>(const triangle<T,2>& triangle, const point2d<T>& point);\
      template inline void closest_point_on_triangle_from_point<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& px, const T& py, const T& pz, T& nx, T& ny, T& nz);\
      template inline point3d<T> closest_point_on_triangle_from_point<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& px, const T& py, const T& pz);\
      template inline point3d<T> closest_point_on_triangle_from_point<T>(const triangle<T,3>& triangle, const T& px, const T& py, const T& pz);\
      template inline point3d<T> closest_point_on_triangle_from_point<T>(const triangle<T,3>& triangle, const point3d<T>& point);\
      template inline void closest_point_on_rectangle_from_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& px, const T& py, T& nx, T& ny);\
      template inline point2d<T> closest_point_on_rectangle_from_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& px, const T& py);\
      template inline point2d<T> closest_point_on_rectangle_from_point<T>(const rectangle<T>& rectangle, const T& px, const T& py);\
      template inline point2d<T> closest_point_on_rectangle_from_point<T>(const rectangle<T>& rectangle, const point2d<T>& point);\
      template inline void closest_point_on_box_from_point<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& px, const T& py, const T& pz, T& nx, T& ny, T& nz);\
      template inline point3d<T> closest_point_on_box_from_point<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& px, const T& py, const T& pz);\
      template inline point3d<T> closest_point_on_box_from_point<T>(const box<T,3>& box, const T& px, const T& py, const T& pz);\
      template inline point3d<T> closest_point_on_box_from_point<T>(const box<T,3>& box, const point3d<T>& point);\
      template inline void closest_point_on_quadix_from_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, const T& px, const T& py, T& nx, T& ny);\
      template inline point2d<T> closest_point_on_quadix_from_point<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, const T& px, const T& py);\
      template inline point2d<T> closest_point_on_quadix_from_point<T>(const quadix<T,2>& quadix, const point2d<T>& point);\
      template inline point2d<T> closest_point_on_circle_from_point<T>(const circle<T>& circle, const point2d<T>& point);\
      template inline point3d<T> closest_point_on_sphere_from_point<T>(const sphere<T>& sphere, const point3d<T>& point);\
      template inline point2d<T> closest_point_on_aabbb_from_point<T>(const rectangle<T>& rectangle, const point2d<T>&point);\
      template inline point2d<T> closest_point_on_circle_from_segment<T>(const circle<T>& circle, const segment<T,2>& segment);\
      template inline point3d<T> closest_point_on_sphere_from_segment<T>(const sphere<T>& sphere, const segment<T,3>& segment);\
      template inline point3d<T> closest_point_on_plane_from_point<T>(const plane<T,3>& plane, const point3d<T>& point);\
      template inline point2d<T> closest_point_on_bezier_from_point<T>(const quadratic_bezier<T,2>& bezier, const point2d<T>& point, const std::size_t& steps);\
      template inline point2d<T> closest_point_on_bezier_from_point<T>(const cubic_bezier<T,2>& bezier, const point2d<T>& point, const std::size_t& steps);\
      template inline point3d<T> closest_point_on_bezier_from_point<T>(const quadratic_bezier<T,3>& bezier, const point3d<T>& point, const std::size_t& steps);\
      template inline point3d<T> closest_point_on_bezier_from_point<T>(const cubic_bezier<T,3>& bezier, const point3d<T>& point, const std::size_t& steps);\
      template inline point2d<T> closest_point_on_circle_from_circle<T>(const circle<T>& circle1, const circle<T>& circle2);\
      template inline point3d<T> closest_point_on_sphere_from_sphere<T>(const sphere<T>& sphere1, const sphere<T>& sphere2);\
      template inline point2d<T> closest_point_on_polygon_from_point<T>(const polygon<T,2>& polygon, const point2d<T>& point);\
      template inline T minimum_distance_from_point_to_segment<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline T minimum_distance_from_point_to_segment<T>(const T& px, const T& py, const T& pz, const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline T minimum_distance_from_point_to_segment<T>(const point2d<T>& point, const segment<T,2>& segment);\
      template inline T minimum_distance_from_point_to_segment<T>(const point3d<T>& point, const segment<T,3>& segment);\
      template inline T minimum_distance_from_point_to_line<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline T minimum_distance_from_point_to_line<T>(const T& px, const T& py, const T& pz, const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline T minimum_distance_from_point_to_line<T>(const point2d<T>& point, const line<T,2>& line);\
      template inline T minimum_distance_from_point_to_line<T>(const point3d<T>& point, const line<T,3>& line);\
      template inline T minimum_distance_from_point_to_triangle<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline T minimum_distance_from_point_to_triangle<T>(const point2d<T>& point, const triangle<T,2>& triangle);\
      template inline T minimum_distance_from_point_to_rectangle<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline T minimum_distance_from_point_to_rectangle<T>(const point2d<T>& point, const rectangle<T>& rectangle);\
      template inline void segment_mid_point<T>(const T&x1, const T&y1, const T&x2, const T&y2, T& midx, T& midy);\
      template inline void segment_mid_point<T>(const segment<T,2>& segment, T& midx, T& midy);\
      template inline point2d<T> segment_mid_point<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline point2d<T> segment_mid_point<T>(const segment<T,2>& segment);\
      template inline void segment_mid_point<T>(const T&x1, const T&y1, const T&z1, const T&x2, const T&y2, const T&z2, T& midx, T& midy, T& midz);\
      template inline void segment_mid_point<T>(const segment<T,3>& segment, T& midx, T& midy, T& midz);\
      template inline point3d<T> segment_mid_point<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline point3d<T> segment_mid_point<T>(const segment<T,3>& segment);\
      template inline void centroid<T>(const T& x1, const T& y1, const T& x2, const T& y2, T& x, T& y);\
      template inline void centroid<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, T& x, T& y, T& z);\
      template inline point2d<T> centroid<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline point2d<T> centroid<T>(const segment<T,2>& segment);\
      template inline void centroid<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, T& x, T& y);\
      template inline void centroid<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, T& x, T& y);\
      template inline void centroid<T>(const triangle<T,2>& triangle, T& x, T& y);\
      template inline void centroid<T>(const triangle<T,3>& triangle, T& x, T& y,T& z);\
      template inline void centroid<T>(const quadix<T,2>& quadix, T& x, T& y);\
      template inline void centroid<T>(const rectangle<T>& rectangle, T& x, T& y);\
      template inline void centroid<T>(const box<T,3>& box, T& x, T& y, T& z);\
      template inline void centroid<T>(const polygon<T,2>& polygon, T& x, T& y);\
      template inline point2d<T> centroid<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline point2d<T> centroid<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4);\
      template inline point2d<T> centroid<T>(const triangle<T,2>& triangle);\
      template inline point3d<T> centroid<T>(const triangle<T,3>& triangle);\
      template inline point2d<T> centroid<T>(const quadix<T,2>& quadix);\
      template inline point2d<T> centroid<T>(const rectangle<T>& rectangle);\
      template inline point3d<T> centroid<T>(const box<T,3>& box);\
      template inline point2d<T> centroid<T>(const polygon<T,2>& polygon);\
      template inline bool point_in_convex_polygon<T>(const T& px, const T& py, const polygon<T,2>& polygon);\
      template inline bool point_in_convex_polygon<T>(const point2d<T>& point, const polygon<T,2>& polygon);\
      template inline bool point_on_polygon_edge<T>(const T& px, const T& py, const polygon<T,2>& polygon);\
      template inline bool point_on_polygon_edge<T>(const point2d<T>& point, const polygon<T,2>& polygon);\
      template inline bool point_in_polygon<T>(const T& px, const T& py, const polygon<T,2>& polygon);\
      template inline bool point_in_polygon<T>(const point2d<T>& point, const polygon<T,2>& polygon);\
      template inline bool point_in_polygon_winding_number<T>(const T& px, const T& py, const polygon<T,2>& polygon);\
      template inline bool point_in_polygon_winding_number<T>(const point2d<T>& point, const polygon<T,2>& polygon);\
      template inline bool convex_quadix<T>(const quadix<T,2>& quadix);\
      template inline bool convex_quadix<T>(const quadix<T,3>& quadix);\
      template inline bool is_convex_polygon<T>(const polygon<T,2>& polygon);\
      template inline polygon<T,2> remove_consecutive_collinear_points<T>(const polygon<T,2>& polygon);\
      template inline void remove_consecutive_collinear_points<T,InputIterator2d,OutputIterator2d>(const InputIterator2d begin, const InputIterator2d end, OutputIterator2d out);\
      template inline bool convex_vertex<T>(const std::size_t& index, const polygon<T,2>& polygon, const int& polygon_orientation);\
      template inline bool collinear_vertex<T>(const std::size_t& index, const polygon<T,2>& polygon);\
      template inline bool vertex_is_ear<T>(const std::size_t& index, const polygon<T,2>& polygon);\
      template inline triangle<T,2> vertex_triangle<T>(const std::size_t& index, const polygon<T,2> polygon);\
      template inline int polygon_orientation<T>(const polygon<T,2>& polygon);\
      template inline bool is_equilateral_triangle<T>(const triangle<T,2>& triangle);\
      template inline bool is_equilateral_triangle<T>(const triangle<T,3>& triangle);\
      template inline bool is_isosceles_triangle<T>(const triangle<T,2>& triangle);\
      template inline bool is_isosceles_triangle<T>(const triangle<T,3>& triangle);\
      template inline bool is_right_triangle<T>(const wykobi::triangle<T,2>& triangle);\
      template inline bool is_right_triangle<T>(const wykobi::triangle<T,3>& triangle);\
      template inline bool are_perspective_triangles<T>(const triangle<T,2>& triangle1, const triangle<T,2>& triangle2);\
      template inline bool are_perspective_triangles<T>(const triangle<T,3>& triangle1, const triangle<T,3>& triangle2);\
      template inline line<T,2> perspectrix<T>(const triangle<T,2>& triangle1, const triangle<T,2>& triangle2);\
      template inline line<T,3> perspectrix<T>(const triangle<T,3>& triangle1, const triangle<T,3>& triangle2);\
      template inline void mirror<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2, T& nx, T& ny);\
      template inline void mirror<T>(const T& px, const T& py, const T& pz, const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, T& nx, T& ny, T& nz);\
      template inline point2d<T> mirror<T>(const point2d<T>& point, const line<T,2>& mirror_axis);\
      template inline segment<T,2> mirror<T>(const segment<T,2>& segment, const line<T,2>& mirror_axis);\
      template inline line<T,2> mirror<T>(const line<T,2>& line, const wykobi::line<T,2>& mirror_axis);\
      template inline rectangle<T> mirror<T>(const rectangle<T>& rectangle, const line<T,2>& mirror_axis);\
      template inline triangle<T,2> mirror<T>(const triangle<T,2>& triangle, const line<T,2>& mirror_axis);\
      template inline quadix<T,2> mirror<T>(const quadix<T,2>& quadix, const line<T,2>& mirror_axis);\
      template inline circle<T> mirror<T>(const circle<T>& circle, const line<T,2>& mirror_axis);\
      template inline polygon<T,2> mirror<T>(const polygon<T,2>& polygon, const line<T,2>& mirror_axis);\
      template inline point3d<T> mirror<T>(const point3d<T>& point, const line<T,3>& mirror_line);\
      template inline segment<T,3> mirror<T>(const segment<T,3>& segment, const line<T,3>& mirror_axis);\
      template inline line<T,3> mirror<T>(const line<T,3>& line, const wykobi::line<T,3>& mirror_axis);\
      template inline box<T,3> mirror<T>(const box<T,3>& box, const line<T,3>& mirror_axis);\
      template inline triangle<T,3> mirror<T>(const triangle<T,3>& triangle, const line<T,3>& mirror_axis);\
      template inline quadix<T,3> mirror<T>(const quadix<T,3>& quadix, const line<T,3>& mirror_axis);\
      template inline sphere<T> mirror<T>(const sphere<T>& sphere, const line<T,3>& mirror_axis);\
      template inline polygon<T,3> mirror<T>(const polygon<T,3>& polygon, const line<T,3>& mirror_axis);\
      template inline point3d<T> mirror<T>(const point3d<T>& point, const plane<T,3>& plane);\
      template inline segment<T,3> mirror<T>(const segment<T,3>& segment, const plane<T,3>& mirror_plane);\
      template inline line<T,3> mirror<T>(const line<T,3>& line, const plane<T,3>& mirror_plane);\
      template inline box<T,3> mirror<T>(const box<T,3>& box, const plane<T,3>& mirror_plane);\
      template inline triangle<T,3> mirror<T>(const triangle<T,3>& triangle, const plane<T,3>& mirror_plane);\
      template inline quadix<T,3> mirror<T>(const quadix<T,3>& quadix, const plane<T,3>& mirror_plane);\
      template inline sphere<T> mirror<T>(const sphere<T>& sphere, const plane<T,3>& mirror_plane);\
      template inline polygon<T,3> mirror<T>(const polygon<T,3>& polygon, const plane<T,3>& mirror_plane);\
      template inline void nonsymmetric_mirror<T>(const T& px, const T& py, const T& x1, const T& y1, const T& x2, const T& y2, const T& ratio, T& nx, T& ny);\
      template inline point2d<T> nonsymmetric_mirror<T>(const point2d<T>& point, const T& ratio, const line<T,2>& line);\
      template inline segment<T,2> nonsymmetric_mirror<T>(const segment<T,2>& segment, const T& ratio, const line<T,2>& line);\
      template inline rectangle<T> nonsymmetric_mirror<T>(const rectangle<T>& rectangle, const T& ratio, const line<T,2>& line);\
      template inline triangle<T,2> nonsymmetric_mirror<T>(const triangle<T,2>& triangle, const T& ratio, const line<T,2>& line);\
      template inline quadix<T,2> nonsymmetric_mirror<T>(const quadix<T,2>& quadix, const T& ratio, const line<T,2>& line);\
      template inline circle<T> nonsymmetric_mirror<T>(const circle<T>& circle, const T& ratio, const line<T,2>& line);\
      template inline polygon<T,2> nonsymmetric_mirror<T>(const polygon<T,2>& polygon, const T& ratio, const line<T,2>& line);\
      template inline point3d<T> nonsymmetric_mirror<T>(const point3d<T>& point, const T& ratio, const plane<T,3>& plane);\
      template inline segment<T,3> nonsymmetric_mirror<T>(const segment<T,3>& segment, const T& ratio, const plane<T,3>& plane);\
      template inline box<T,3> nonsymmetric_mirror<T>(const box<T,3>& box, const T& ratio, const plane<T,3>& plane);\
      template inline triangle<T,3> nonsymmetric_mirror<T>(const triangle<T,3>& triangle, const T& ratio, const plane<T,3>& plane);\
      template inline quadix<T,3> nonsymmetric_mirror<T>(const quadix<T,3>& quadix, const T& ratio, const plane<T,3>& plane);\
      template inline sphere<T> nonsymmetric_mirror<T>(const sphere<T>& sphere, const T& ratio, const plane<T,3>& plane);\
      template inline polygon<T,3> nonsymmetric_mirror<T>(const polygon<T,3>& polygon, const T& ratio, const plane<T,3>& plane);\
      template inline point2d<T> invert_point<T>(const point2d<T>& point, const circle<T>& circle);\
      template inline point3d<T> invert_point<T>(const point3d<T>& point, const sphere<T>& sphere);\
      template inline point2d<T> antipodal_point<T>(const point2d<T>& point, const circle<T>& circle);\
      template inline point3d<T> antipodal_point<T>(const point3d<T>& point, const sphere<T>& sphere);\
      template inline T distance<T>(const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline T distance<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline T distance<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline T distance<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline T distance<T>(const curve_point<T,2>& point1, const curve_point<T,2>& point2);\
      template inline T distance<T>(const curve_point<T,3>& point1, const curve_point<T,3>& point2);\
      template inline T distance<T>(const point2d<T>& point, const segment<T,2>& segment);\
      template inline T distance<T>(const point3d<T>& point, const segment<T,3>& segment);\
      template inline T distance<T>(const point2d<T>& point, const rectangle<T>& rectangle);\
      template inline T distance<T>(const point2d<T>& point, const triangle<T,2>& triangle);\
      template inline T distance<T>(const point2d<T>& point, const quadix<T,2>& quadix);\
      template inline T distance<T>(const point2d<T>& point, const ray<T,2>& ray);\
      template inline T distance<T>(const point3d<T>& point, const ray<T,3>& ray);\
      template inline T distance<T>(const point3d<T>& point, const plane<T,3>& plane);\
      template inline T distance<T>(const line<T,2>& line1, const line<T,2>& line2);\
      template inline T distance<T>(const line<T,3>& line1, const line<T,3>& line2);\
      template inline T distance<T>(const segment<T,2>& segment1, const segment<T,2>& segment2);\
      template inline T distance<T>(const segment<T,3>& segment1, const segment<T,3>& segment2);\
      template inline T distance<T>(const segment<T,2>& segment);\
      template inline T distance<T>(const segment<T,3>& segment);\
      template inline T distance<T>(const segment<T,2>& segment, const triangle<T,2>& triangle);\
      template inline T distance<T>(const segment<T,3>& segment, const triangle<T,3>& triangle);\
      template inline T distance<T>(const segment<T,2>& segment, const rectangle<T>& rectangle);\
      template inline T distance<T>(const segment<T,2>& segment, const circle<T>& circle);\
      template inline T distance<T>(const triangle<T,2>& triangle1, const triangle<T,2>& triangle2);\
      template inline T distance<T>(const triangle<T,2>& triangle, const rectangle<T>& rectangle);\
      template inline T distance<T>(const rectangle<T>& rectangle1, const rectangle<T>& rectangle2);\
      template inline T distance<T>(const triangle<T,2>& triangle, const circle<T>& circle);\
      template inline T distance<T>(const rectangle<T>& rectangle, const circle<T>& circle);\
      template inline T distance<T>(const point2d<T>& point, const circle<T>& circle);\
      template inline T distance<T>(const circle<T>& circle1, const circle<T>& circle2);\
      template inline T distance<T>(const sphere<T>& sphere1, const sphere<T>& sphere2);\
      template inline T lay_distance<T>(const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline T lay_distance<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline T lay_distance<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline T lay_distance<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline T lay_distance<T>(const point2d<T>& point, const triangle<T,2>& triangle);\
      template inline T lay_distance<T>(const point2d<T>& point, const quadix<T,2>& triangle);\
      template inline T lay_distance<T>(const point2d<T>& point, const ray<T,2>& ray);\
      template inline T lay_distance<T>(const point3d<T>& point, const ray<T,3>& ray);\
      template inline T lay_distance<T>(const point3d<T>& point, const plane<T,3>& plane);\
      template inline T lay_distance<T>(const segment<T,2>& segment1, const segment<T,2>& segment2);\
      template inline T lay_distance<T>(const segment<T,3>& segment1, const segment<T,3>& segment2);\
      template inline T lay_distance<T>(const line<T,3>& line1, const line<T,3>& line2);\
      template inline T lay_distance<T>(const segment<T,2>& segment);\
      template inline T lay_distance<T>(const segment<T,3>& segment);\
      template inline T lay_distance<T>(const segment<T,2>& segment, const triangle<T,2>& triangle);\
      template inline T lay_distance<T>(const segment<T,3>& segment, const triangle<T,3>& triangle);\
      template inline T manhattan_distance<T>(const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline T manhattan_distance<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline T manhattan_distance<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline T manhattan_distance<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline T manhattan_distance<T>(const point2d<T>& point, const ray<T,2>& ray);\
      template inline T manhattan_distance<T>(const point3d<T>& point, const ray<T,3>& ray);\
      template inline T manhattan_distance<T>(const segment<T,2>& segment);\
      template inline T manhattan_distance<T>(const segment<T,3>& segment);\
      template inline T manhattan_distance<T>(const circle<T>& circle1, const circle<T>& circle2);\
      template inline T chebyshev_distance<T>(const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline T chebyshev_distance<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline T chebyshev_distance<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline T chebyshev_distance<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline T chebyshev_distance<T>(const segment<T,2>& segment);\
      template inline T chebyshev_distance<T>(const segment<T,3>& segment);\
      template inline T chebyshev_distance<T>(const circle<T>& circle1, const circle<T>& circle2);\
      template inline T inverse_chebyshev_distance<T>(const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline T inverse_chebyshev_distance<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline T inverse_chebyshev_distance<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline T inverse_chebyshev_distance<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline T inverse_chebyshev_distance<T>(const segment<T,2>& segment);\
      template inline T inverse_chebyshev_distance<T>(const segment<T,3>& segment);\
      template inline T inverse_chebyshev_distance<T>(const circle<T>& circle1, const circle<T>& circle2);\
      template inline point2d<T> minkowski_sum<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline polygon<T,2> minkowski_sum<T>(const rectangle<T>& rectangle1, const rectangle<T>& rectangle2);\
      template inline polygon<T,2> minkowski_sum<T>(const triangle<T,2>& triangle1, const triangle<T,2>& triangle2);\
      template inline polygon<T,2> minkowski_sum<T>(const quadix<T,2>& quadix1, const quadix<T,2>& quadix2);\
      template inline polygon<T,2> minkowski_sum<T>(const circle<T>& triangle, const circle<T>& circle);\
      template inline polygon<T,2> minkowski_sum<T>(const triangle<T,2>& triangle, const rectangle<T>& rectangle);\
      template inline polygon<T,2> minkowski_sum<T>(const triangle<T,2>& triangle, const quadix<T,2>& quadix);\
      template inline polygon<T,2> minkowski_sum<T>(const triangle<T,2>& triangle, const circle<T>& circle);\
      template inline polygon<T,2> minkowski_sum<T>(const quadix<T,2>& quadix, const circle<T>& circle);\
      template inline polygon<T,2> minkowski_sum<T>(const quadix<T,2>& quadix, const rectangle<T>& rectangle);\
      template inline polygon<T,2> minkowski_sum<T>(const rectangle<T>& rectangle, const circle<T>& circle);\
      template inline polygon<T,2> minkowski_sum<T>(const polygon<T,2>& polygon1, const polygon<T,2>& polygon2);\
      template inline point2d<T> minkowski_difference<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline polygon<T,2> minkowski_difference<T>(const rectangle<T>& rectangle1, const rectangle<T>& rectangle2);\
      template inline polygon<T,2> minkowski_difference<T>(const triangle<T,2>& triangle1, const triangle<T,2>& triangle2);\
      template inline polygon<T,2> minkowski_difference<T>(const quadix<T,2>& quadix1, const quadix<T,2>& quadix2);\
      template inline polygon<T,2> minkowski_difference<T>(const circle<T>& triangle, const circle<T>& circle);\
      template inline polygon<T,2> minkowski_difference<T>(const triangle<T,2>& triangle, const rectangle<T>& rectangle);\
      template inline polygon<T,2> minkowski_difference<T>(const triangle<T,2>& triangle, const quadix<T,2>& quadix);\
      template inline polygon<T,2> minkowski_difference<T>(const triangle<T,2>& triangle, const circle<T>& circle);\
      template inline polygon<T,2> minkowski_difference<T>(const quadix<T,2>& quadix, const circle<T>& circle);\
      template inline polygon<T,2> minkowski_difference<T>(const quadix<T,2>& quadix, const rectangle<T>& rectangle);\
      template inline polygon<T,2> minkowski_difference<T>(const rectangle<T>& rectangle, const circle<T>& circle);\
      template inline polygon<T,2> minkowski_difference<T>(const polygon<T,2>& polygon1, const polygon<T,2>& polygon2);\
      template inline T distance_segment_to_segment<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline T distance_segment_to_segment<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4);\
      template inline T lay_distance_segment_to_segment<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline T lay_distance_segment_to_segment<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4);\
      template inline T distance_line_to_line<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline T distance_line_to_line<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4);\
      template inline T lay_distance_line_to_line<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline T lay_distance_line_to_line<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4);\
      template inline T lay_distance_from_point_to_circle_center<T>(const point2d<T>& point, const circle<T>& circle);\
      template inline T lay_distance_from_point_to_sphere_center<T>(const point3d<T>& point, const sphere<T>& sphere);\
      template inline T distance_from_point_to_circle_center<T>(const point2d<T>& point, const circle<T>& circle);\
      template inline T distance_from_point_to_sphere_center<T>(const point3d<T>& point, const sphere<T>& sphere);\
      template inline T span_length(const rectangle<T>& rect);\
      template inline T span_length(const box<T,3>& box);\
      template inline void project_point_t<T>(const T& srcx, const T& srcy, const T& destx, const T& desty, const T& t, T& nx, T& ny);\
      template inline void project_point_t<T>(const T& srcx, const T& srcy, const T& srcz, const T& destx, const T& desty, const T& destz, const T& t, T& nx, T& ny, T& nz);\
      template inline void project_point<T>(const T& srcx, const T& srcy, const T& destx, const T& desty, const T& dist, T& nx, T& ny);\
      template inline void project_point<T>(const T& srcx, const T& srcy, const T& srcz, const T& destx, const T& desty, const T& destz, const T& dist, T& nx, T& ny, T& nz);\
      template inline void project_point<T>(const T& px, const T& py, const T& angle, const T& distance, T& nx, T& ny);\
      template inline void project_point0 <T>(const T& px, const T& py, const T& distance, T& nx, T& ny);\
      template inline void project_point45 <T>(const T& px, const T& py, const T& distance, T& nx, T& ny);\
      template inline void project_point90 <T>(const T& px, const T& py, const T& distance, T& nx, T& ny);\
      template inline void project_point135<T>(const T& px, const T& py, const T& distance, T& nx, T& ny);\
      template inline void project_point180<T>(const T& px, const T& py, const T& distance, T& nx, T& ny);\
      template inline void project_point225<T>(const T& px, const T& py, const T& distance, T& nx, T& ny);\
      template inline void project_point270<T>(const T& px, const T& py, const T& distance, T& nx, T& ny);\
      template inline void project_point315<T>(const T& px, const T& py, const T& distance, T& nx, T& ny);\
      template inline point2d<T> project_point_t<T>(const point2d<T>& source_point, const point2d<T>& destination_point, const T& t);\
      template inline point3d<T> project_point_t<T>(const point3d<T>& source_point, const point3d<T>& destination_point, const T& t);\
      template inline point2d<T> project_point<T>(const point2d<T>& source_point, const point2d<T>& destination_point, const T& distance);\
      template inline point3d<T> project_point<T>(const point3d<T>& source_point, const point3d<T>& destination_point, const T& distance);\
      template inline point2d<T> project_point<T>(const point2d<T>& point, const T& angle, const T& distance);\
      template inline point2d<T> project_point0 <T>(const point2d<T>& point, const T& distance);\
      template inline point2d<T> project_point45 <T>(const point2d<T>& point, const T& distance);\
      template inline point2d<T> project_point90 <T>(const point2d<T>& point, const T& distance);\
      template inline point2d<T> project_point135<T>(const point2d<T>& point, const T& distance);\
      template inline point2d<T> project_point180<T>(const point2d<T>& point, const T& distance);\
      template inline point2d<T> project_point225<T>(const point2d<T>& point, const T& distance);\
      template inline point2d<T> project_point270<T>(const point2d<T>& point, const T& distance);\
      template inline point2d<T> project_point315<T>(const point2d<T>& point, const T& distance);\
      template inline point2d<T> project_object<T>(const point2d<T>& point, const T& angle, const T& distance);\
      template inline segment<T,2> project_object<T>(const segment<T,2>& segment, const T& angle, const T& distance);\
      template inline triangle<T,2> project_object<T>(const triangle<T,2>& triangle, const T& angle, const T& distance);\
      template inline quadix<T,2> project_object<T>(const quadix<T,2>& quadix, const T& angle, const T& distance);\
      template inline circle<T> project_object<T>(const circle<T>& circle, const T& angle, const T& distance);\
      template inline polygon<T,2> project_object<T>(const polygon<T,2>& polygon, const T& angle, const T& distance);\
      template inline segment<T,2> project_onto_axis<T>(const point2d<T>& point, const line<T,2>& axis);\
      template inline segment<T,2> project_onto_axis<T>(const triangle<T,2>& triangle, const line<T,2>& axis);\
      template inline segment<T,2> project_onto_axis<T>(const rectangle<T>& rectangle, const line<T,2>& axis);\
      template inline segment<T,2> project_onto_axis<T>(const quadix<T,2>& quadix, const line<T,2>& axis);\
      template inline segment<T,2> project_onto_axis<T>(const circle<T>& circle, const line<T,2>& axis);\
      template inline segment<T,2> project_onto_axis<T>(const polygon<T,2>& polygon, const line<T,2>& axis);\
      template inline segment<T,3> project_onto_axis<T>(const point3d<T>& point, const line<T,3>& axis);\
      template inline segment<T,3> project_onto_axis<T>(const triangle<T,3>& triangle, const line<T,3>& axis);\
      template inline segment<T,3> project_onto_axis<T>(const box<T,3>& box, const line<T,3>& axis);\
      template inline segment<T,3> project_onto_axis<T>(const quadix<T,3>& quadix, const line<T,3>& axis);\
      template inline segment<T,3> project_onto_axis<T>(const sphere<T>& sphere, const line<T,3>& axis);\
      template inline segment<T,3> project_onto_axis<T>(const polygon<T,3>& polygon, const line<T,3>& axis);\
      template inline void calculate_bezier_coefficients<T>(const quadratic_bezier<T,2>& bezier, T& ax, T& bx, T& ay, T& by);\
      template inline void calculate_bezier_coefficients<T>(const quadratic_bezier<T,3>& bezier, T& ax, T& bx, T& ay, T& by, T& az, T& bz);\
      template inline void calculate_bezier_coefficients<T>(const cubic_bezier<T,2>& bezier, T& ax, T& bx, T& cx, T& ay, T& by, T& cy);\
      template inline void calculate_bezier_coefficients<T>(const cubic_bezier<T,3>& bezier, T& ax, T& bx, T& cx, T& ay, T& by, T& cy, T& az, T& bz, T& cz);\
      template inline void calculate_bezier_coefficients<T>(const quadratic_bezier<T,2>& bezier, bezier_coefficients<T,2,eQuadraticBezier>& coeffs);\
      template inline void calculate_bezier_coefficients<T>(const quadratic_bezier<T,3>& bezier, bezier_coefficients<T,3,eQuadraticBezier>& coeffs);\
      template inline void calculate_bezier_coefficients<T>(const cubic_bezier<T,2>& bezier, bezier_coefficients<T,2,eCubicBezier>& coeffs);\
      template inline void calculate_bezier_coefficients<T>(const cubic_bezier<T,3>& bezier, bezier_coefficients<T,3,eCubicBezier>& coeffs);\
      template inline point2d<T> create_point_on_bezier<T>(const point2d<T>& start_point, const T& ax, const T& bx, const T& ay, const T& by, const T& t);\
      template inline point3d<T> create_point_on_bezier<T>(const point3d<T>& start_point, const T& ax, const T& bx, const T& ay, const T& by, const T& az, const T& bz, const T& t);\
      template inline point2d<T> create_point_on_bezier<T>(const point2d<T>& start_point, const T& ax, const T& bx, const T& cx, const T& ay, const T& by, const T& cy, const T& t);\
      template inline point3d<T> create_point_on_bezier<T>(const point3d<T>& start_point, const T& ax, const T& bx, const T& cx, const T& ay, const T& by, const T& cy, const T& az, const T& bz, const T& cz, const T& t);\
      template inline point2d<T> create_point_on_bezier<T>(const point2d<T>& start_point, const bezier_coefficients<T,2,eQuadraticBezier>& coeffs, const T& t);\
      template inline point3d<T> create_point_on_bezier<T>(const point3d<T>& start_point, const bezier_coefficients<T,3,eQuadraticBezier>& coeffs, const T& t);\
      template inline point2d<T> create_point_on_bezier<T>(const point2d<T>& start_point, const bezier_coefficients<T,2,eCubicBezier>& coeffs, const T& t);\
      template inline point3d<T> create_point_on_bezier<T>(const point3d<T>& start_point, const bezier_coefficients<T,3,eCubicBezier>& coeffs, const T& t);\
      template inline void generate_bezier<T,OutputIterator2d>(const quadratic_bezier<T,2>& bezier, OutputIterator2d out, const std::size_t& point_count);\
      template inline void generate_bezier<T,OutputIterator3d>(const quadratic_bezier<T,3>& bezier, OutputIterator3d out, const std::size_t& point_count);\
      template inline void generate_bezier<T,OutputIterator2d>(const cubic_bezier<T,2>& bezier, OutputIterator2d out, const std::size_t& point_count);\
      template inline void generate_bezier<T,OutputIterator3d>(const cubic_bezier<T,3>& bezier, OutputIterator3d out, const std::size_t& point_count);\
      template inline T bezier_curve_length<T>(const quadratic_bezier<T,2>& bezier, const std::size_t& point_count);\
      template inline T bezier_curve_length<T>(const quadratic_bezier<T,3>& bezier, const std::size_t& point_count);\
      template inline T bezier_curve_length<T>(const cubic_bezier<T,2>& bezier, const std::size_t& point_count);\
      template inline T bezier_curve_length<T>(const cubic_bezier<T,3>& bezier, const std::size_t& point_count);\
      template inline triangle<T,2> bezier_convex_hull<T>(const quadratic_bezier<T,2>& bezier);\
      template inline quadix<T,2> bezier_convex_hull<T>(const cubic_bezier<T,2>& bezier);\
      template inline segment<T,2> center_at_location<T>(const segment<T,2>& segment, const T& x, const T& y);\
      template inline segment<T,3> center_at_location<T>(const segment<T,3>& segment, const T& x, const T& y, const T& z);\
      template inline triangle<T,2> center_at_location<T>(const triangle<T,2>& triangle, const T& x, const T& y);\
      template inline rectangle<T> center_at_location<T>(const rectangle<T>& rectangle, const T& x, const T& y);\
      template inline box<T,3> center_at_location<T>(const box<T,3>& box, const T& x, const T& y, const T& z);\
      template inline quadix<T,2> center_at_location<T>(const quadix<T,2>& quadix, const T& x, const T& y);\
      template inline circle<T> center_at_location<T>(const circle<T>& circle, const T& x, const T& y);\
      template inline polygon<T,2> center_at_location<T>(const polygon<T,2>& polygon, const T& x, const T& y);\
      template inline segment<T,2> center_at_location<T>(const segment<T,2>& segment, const point2d<T>& center_point);\
      template inline segment<T,3> center_at_location<T>(const segment<T,3>& segment, const point3d<T>& center_point);\
      template inline triangle<T,2> center_at_location<T>(const triangle<T,2>& triangle, const point2d<T>& center_point);\
      template inline rectangle<T> center_at_location<T>(const rectangle<T>& rectangle, const point2d<T>& center_point);\
      template inline box<T,3> center_at_location<T>(const box<T,3>& box, const point3d<T>& center_point);\
      template inline quadix<T,2> center_at_location<T>(const quadix<T,2>& quadix, const point2d<T>& center_point);\
      template inline circle<T> center_at_location<T>(const circle<T>& circle, const point2d<T>& center_point);\
      template inline polygon<T,2> center_at_location<T>(const polygon<T,2>& polygon, const point2d<T>& center_point);\
      template inline void shorten_segment<T>(T& x1, T& y1, T& x2, T& y2, const T& amount);\
      template inline void shorten_segment<T>(T& x1, T& y1, T& z1, T& x2, T& y2, T& z2, const T& amount);\
      template inline segment<T,2> shorten_segment<T>(const segment<T,2>& segment, const T& amount);\
      template inline segment<T,3> shorten_segment<T>(const segment<T,3>& segment, const T& amount);\
      template inline void lengthen_segment<T>(T& x1, T& y1, T& x2, T& y2, const T& amount);\
      template inline void lengthen_segment<T>(T& x1, T& y1, T& z1, T& x2, T& y2, T& z2, const T& amount);\
      template inline segment<T,2> lengthen_segment<T>(const segment<T,2>& segment, const T& amount);\
      template inline segment<T,3> lengthen_segment<T>(const segment<T,3>& segment, const T& amount);\
      template inline int out_code<T>(const point2d<T>& point, const rectangle<T>& rectangle);\
      template inline bool clip<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4, T& cx1,T& cy1, T& cx2,T& cy2);\
      template inline bool clip<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4, T& cx1,T& cy1,T& cz1, T& cx2,T& cy2,T& cz2);\
      template inline bool clip<T>(const segment<T,2>& src_segment, const rectangle<T>& rectangle, segment<T,2>& csegment);\
      template inline bool clip<T>(const segment<T,2>& src_segment, const triangle<T,2>& triangle,segment<T,2>& csegment);\
      template inline bool clip<T>(const segment<T,2>& src_segment, const quadix<T,2>&quadix, segment<T,2>& csegment);\
      template inline bool clip<T>(const segment<T,2>& src_segment, const circle<T>& circle, segment<T,2>& csegment);\
      template inline bool clip<T>(const rectangle<T>& rectangle1, const rectangle<T>& rectangle2, rectangle<T>& crectangle);\
      template inline bool clip<T>(const box<T,3>& box1, const box<T,3>& box2, box<T,3>& cbox);\
      template inline T area<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline T area<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3);\
      template inline T area<T>(const triangle<T,2>& triangle);\
      template inline T area<T>(const triangle<T,3>& triangle);\
      template inline T area<T>(const quadix<T,2>& quadix);\
      template inline T area<T>(const quadix<T,3>& quadix);\
      template inline T area<T>(const rectangle<T>& rectangle);\
      template inline T area<T>(const circle<T>& circle);\
      template inline T area<T>(const polygon<T,2>& polygon);\
      template inline T perimeter<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline T perimeter<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3);\
      template inline T perimeter<T>(const triangle<T,2>& triangle);\
      template inline T perimeter<T>(const triangle<T,3>& triangle);\
      template inline T perimeter<T>(const quadix<T,2>& quadix);\
      template inline T perimeter<T>(const quadix<T,3>& quadix);\
      template inline T perimeter<T>(const rectangle<T>& rectangle);\
      template inline T perimeter<T>(const circle<T>& circle);\
      template inline T perimeter<T>(const polygon<T,2>& polygon);\
      template inline void rotate<T>(const T& rotation_angle, const T& x, const T& y, T& nx, T& ny);\
      template inline void rotate<T>(const T& rotation_angle, const T& x, const T& y, const T& ox, const T& oy, T& nx, T& ny);\
      template inline point2d<T> rotate<T>(const T& rotation_angle, const point2d<T>& point);\
      template inline point2d<T> rotate<T>(const T& rotation_angle, const point2d<T>& point, const point2d<T>& opoint);\
      template inline segment<T,2> rotate<T>(const T& rotation_angle, const segment<T,2>& segment);\
      template inline segment<T,2> rotate<T>(const T& rotation_angle, const segment<T,2>& segment, const point2d<T>& opoint);\
      template inline triangle<T,2> rotate<T>(const T& rotation_angle, const triangle<T,2>& triangle);\
      template inline triangle<T,2> rotate<T>(const T& rotation_angle, const triangle<T,2>& triangle, const point2d<T>& opoint);\
      template inline quadix<T,2> rotate<T>(const T& rotation_angle, const quadix<T,2>& quadix);\
      template inline quadix<T,2> rotate<T>(const T& rotation_angle, const quadix<T,2>& quadix, const point2d<T>& opoint);\
      template inline polygon<T,2> rotate<T>(const T& rotation_angle, const polygon<T,2>& polygon);\
      template inline polygon<T,2> rotate<T>(const T& rotation_angle, const polygon<T,2>& polygon, const point2d<T>& opoint);\
      template inline void rotate<T>(const T& rx, const T& ry, const T& rz, const T& x, const T& y, const T& z, T& nx, T& ny, T& nz);\
      template inline void rotate<T>(const T& rx, const T& ry, const T& rz, const T& x, const T& y, const T& z, const T& ox, const T& oy, const T& oz, T& nx, T& ny, T& nz);\
      template inline point3d<T> rotate<T>(const T& rx, const T& ry, const T& rz, const point3d<T>& point);\
      template inline point3d<T> rotate<T>(const T& rx, const T& ry, const T& rz, const point3d<T>& point, const point3d<T>& opoint);\
      template inline segment<T,3> rotate<T>(const T& rx, const T& ry, const T& rz, const segment<T,3>& segment);\
      template inline segment<T,3> rotate<T>(const T& rx, const T& ry, const T& rz, const segment<T,3>& segment, const point3d<T>& opoint);\
      template inline triangle<T,3> rotate<T>(const T& rx, const T& ry, const T& rz, const triangle<T,3>& triangle);\
      template inline triangle<T,3> rotate<T>(const T& rx, const T& ry, const T& rz, const triangle<T,3>& triangle, const point3d<T>& opoint);\
      template inline quadix<T,3> rotate<T>(const T& rx, const T& ry, const T& rz, const quadix<T,3>& quadix);\
      template inline quadix<T,3> rotate<T>(const T& rx, const T& ry, const T& rz, const quadix<T,3>& quadix, const point3d<T>& opoint);\
      template inline polygon<T,3> rotate<T>(const T& rx, const T& ry, const T& rz, const polygon<T,3>& polygon);\
      template inline polygon<T,3> rotate<T>(const T& rx, const T& ry, const T& rz, const polygon<T,3>& polygon, const point3d<T>& opoint);\
      template inline void fast_rotate<T>(const int rotation_angle, const T& x, const T& y, T& nx, T& ny);\
      template inline void fast_rotate<T>(const int rotation_angle, const T& x, const T& y, const T& ox, const T& oy, T& nx, T& ny);\
      template inline point2d<T> fast_rotate<T>(const int rotation_angle, const point2d<T>& point);\
      template inline point2d<T> fast_rotate<T>(const int rotation_angle, const point2d<T>& point, const point2d<T>& opoint);\
      template inline segment<T,2> fast_rotate<T>(const int rotation_angle, const segment<T,2>& segment);\
      template inline segment<T,2> fast_rotate<T>(const int rotation_angle, const segment<T,2>& segment, const point2d<T>& opoint);\
      template inline triangle<T,2> fast_rotate<T>(const int rotation_angle, const triangle<T,2>& triangle);\
      template inline triangle<T,2> fast_rotate<T>(const int rotation_angle, const triangle<T,2>& triangle, const point2d<T>& opoint);\
      template inline quadix<T,2> fast_rotate<T>(const int rotation_angle, const quadix<T,2>& quadix);\
      template inline quadix<T,2> fast_rotate<T>(const int rotation_angle, const quadix<T,2>& quadix, const point2d<T>& opoint);\
      template inline polygon<T,2> fast_rotate<T>(const int rotation_angle, const polygon<T,2>& polygon);\
      template inline polygon<T,2> fast_rotate<T>(const int rotation_angle, const polygon<T,2>& polygon, const point2d<T>& opoint);\
      template inline void fast_rotate<T>(const int rx, const int ry, const int rz, const T&x, const T&y, const T&z, T& nx, T& ny, T& nz);\
      template inline void fast_rotate<T>(const int rx, const int ry, const int rz, const T&x, const T&y, const T&z, const T& ox, const T& oy, const T& oz, T& nx, T& ny, T& nz);\
      template inline point3d<T> fast_rotate<T>(const int rx, const int ry, const int rz, const point3d<T>& point);\
      template inline point3d<T> fast_rotate<T>(const int rx, const int ry, const int rz, const point3d<T>& point, const point3d<T>& opoint);\
      template inline segment<T,3> fast_rotate<T>(const int rx, const int ry, const int rz, const segment<T,3>& segment);\
      template inline segment<T,3> fast_rotate<T>(const int rx, const int ry, const int rz, const segment<T,3>& segment, const point3d<T>& opoint);\
      template inline triangle<T,3> fast_rotate<T>(const int rx, const int ry, const int rz, const triangle<T,3>& triangle);\
      template inline triangle<T,3> fast_rotate<T>(const int rx, const int ry, const int rz, const triangle<T,3>& triangle, const point3d<T>& opoint);\
      template inline quadix<T,3> fast_rotate<T>(const int rx, const int ry, const int rz, const quadix<T,3>& quadix);\
      template inline quadix<T,3> fast_rotate<T>(const int rx, const int ry, const int rz, const quadix<T,3>& quadix, const point3d<T>& opoint);\
      template inline polygon<T,3> fast_rotate<T>(const int rx, const int ry, const int rz, const polygon<T,3>& polygon);\
      template inline polygon<T,3> fast_rotate<T>(const int rx, const int ry, const int rz, const polygon<T,3>& polygon, const point3d<T>& opoint);\
      template inline point2d<T> translate<T>(const T& dx, const T& dy, const point2d<T>& point);\
      template inline line<T,2> translate<T>(const T& dx, const T& dy, const line<T,2>& line);\
      template inline segment<T,2> translate<T>(const T& dx, const T& dy, const segment<T,2>& segment);\
      template inline triangle<T,2> translate<T>(const T& dx, const T& dy, const triangle<T,2>& triangle);\
      template inline quadix<T,2> translate<T>(const T& dx, const T& dy, const quadix<T,2>& quadix);\
      template inline rectangle<T> translate<T>(const T& dx, const T& dy, const rectangle<T>& rectangle);\
      template inline circle<T> translate<T>(const T& dx, const T& dy, const circle<T>& circle);\
      template inline polygon<T,2> translate<T>(const T& dx, const T& dy, const polygon<T,2>& polygon);\
      template inline point2d<T> translate<T>(const T& delta, const point2d<T>& point);\
      template inline line<T,2> translate<T>(const T& delta, const line<T,2>& line);\
      template inline segment<T,2> translate<T>(const T& delta, const segment<T,2>& segment);\
      template inline triangle<T,2> translate<T>(const T& delta, const triangle<T,2>& triangle);\
      template inline quadix<T,2> translate<T>(const T& delta, const quadix<T,2>& quadix);\
      template inline rectangle<T> translate<T>(const T& delta, const rectangle<T>& rectangle);\
      template inline circle<T> translate<T>(const T& delta, const circle<T>& circle);\
      template inline polygon<T,2> translate<T>(const T& delta, const polygon<T,2>& polygon);\
      template inline point2d<T> translate<T>(const vector2d<T>& v, const point2d<T>& point);\
      template inline line<T,2> translate<T>(const vector2d<T>& v, const line<T,2>& line);\
      template inline segment<T,2> translate<T>(const vector2d<T>& v, const segment<T,2>& segment);\
      template inline triangle<T,2> translate<T>(const vector2d<T>& v, const triangle<T,2>& triangle);\
      template inline quadix<T,2> translate<T>(const vector2d<T>& v, const quadix<T,2>& quadix);\
      template inline rectangle<T> translate<T>(const vector2d<T>& v, const rectangle<T>& rectangle);\
      template inline circle<T> translate<T>(const vector2d<T>& v, const circle<T>& circle);\
      template inline polygon<T,2> translate<T>(const vector2d<T>& v, const polygon<T,2>& polygon);\
      template inline point3d<T> translate<T>(const T& dx, const T& dy, const T& dz, const point3d<T>& point);\
      template inline line<T,3> translate<T>(const T& dx, const T& dy, const T& dz, const line<T,3>& line);\
      template inline segment<T,3> translate<T>(const T& dx, const T& dy, const T& dz, const segment<T,3>& segment);\
      template inline triangle<T,3> translate<T>(const T& dx, const T& dy, const T& dz, const triangle<T,3>& triangle);\
      template inline quadix<T,3> translate<T>(const T& dx, const T& dy, const T& dz, const quadix<T,3>& quadix);\
      template inline box<T,3> translate<T>(const T& dx, const T& dy, const T& dz, const box<T,3>& box);\
      template inline sphere<T> translate<T>(const T& dx, const T& dy, const T& dz, const sphere<T>& sphere);\
      template inline polygon<T,3> translate<T>(const T& dx, const T& dy, const T& dz, const polygon<T,3>& polygon);\
      template inline point3d<T> translate<T>(const T& delta, const point3d<T>& point);\
      template inline line<T,3> translate<T>(const T& delta, const line<T,3>& line);\
      template inline segment<T,3> translate<T>(const T& delta, const segment<T,3>& segment);\
      template inline triangle<T,3> translate<T>(const T& delta, const triangle<T,3>& triangle);\
      template inline quadix<T,3> translate<T>(const T& delta, const quadix<T,3>& quadix);\
      template inline box<T,3> translate<T>(const T& delta, const box<T,3>& box);\
      template inline sphere<T> translate<T>(const T& delta, const sphere<T>& sphere);\
      template inline polygon<T,3> translate<T>(const T& delta, const polygon<T,3>& polygon);\
      template inline point3d<T> translate<T>(const vector3d<T>& v, const point3d<T>& point);\
      template inline line<T,3> translate<T>(const vector3d<T>& v, const line<T,3>& line);\
      template inline segment<T,3> translate<T>(const vector3d<T>& v, const segment<T,3>& segment);\
      template inline triangle<T,3> translate<T>(const vector3d<T>& v, const triangle<T,3>& triangle);\
      template inline quadix<T,3> translate<T>(const vector3d<T>& v, const quadix<T,3>& quadix);\
      template inline box<T,3> translate<T>(const vector3d<T>& v, const box<T,3>& box);\
      template inline sphere<T> translate<T>(const vector3d<T>& v, const sphere<T>& sphere);\
      template inline polygon<T,3> translate<T>(const vector3d<T>& v, const polygon<T,3>& polygon);\
      template inline point2d<T> scale<T>(const T& dx, const T& dy, const point2d<T>& point);\
      template inline line<T,2> scale<T>(const T& dx, const T& dy, const line<T,2>& line);\
      template inline segment<T,2> scale<T>(const T& dx, const T& dy, const segment<T,2>& segment);\
      template inline triangle<T,2> scale<T>(const T& dx, const T& dy, const triangle<T,2>& triangle);\
      template inline quadix<T,2> scale<T>(const T& dx, const T& dy, const quadix<T,2>& quadix);\
      template inline rectangle<T> scale<T>(const T& dx, const T& dy, const rectangle<T>& rectangle);\
      template inline circle<T> scale<T>(const T& dr, const circle<T>& circle);\
      template inline polygon<T,2> scale<T>(const T& dx, const T& dy, const polygon<T,2>& polygon);\
      template inline point3d<T> scale<T>(const T& dx, const T& dy, const T& dz, const point3d<T>& point);\
      template inline line<T,3> scale<T>(const T& dx, const T& dy, const T& dz, const line<T,3>& line);\
      template inline segment<T,3> scale<T>(const T& dx, const T& dy, const T& dz, const segment<T,3>& segment);\
      template inline triangle<T,3> scale<T>(const T& dx, const T& dy, const T& dz, const triangle<T,3>& triangle);\
      template inline quadix<T,3> scale<T>(const T& dx, const T& dy, const T& dz, const quadix<T,3>& quadix);\
      template inline box<T,3> scale<T>(const T& dx, const T& dy, const T& dz, const box<T,3>& box);\
      template inline sphere<T> scale<T>(const T& dr, const sphere<T>& sphere);\
      template inline polygon<T,3> scale<T>(const T& dx, const T& dy, const T& dz, const polygon<T,3>& polygon);\
      template inline rectangle<T> aabb<T>(const segment<T,2>& segment);\
      template inline rectangle<T> aabb<T>(const triangle<T,2>& triangle);\
      template inline rectangle<T> aabb<T>(const rectangle<T>& rectangle);\
      template inline rectangle<T> aabb<T>(const quadix<T,2>& quadix);\
      template inline rectangle<T> aabb<T>(const circle<T>& circle);\
      template inline rectangle<T> aabb<T>(const polygon<T,2>& polygon);\
      template inline void aabb<T>(const segment<T,2>& segment,T& x1, T& y1, T& x2, T& y2);\
      template inline void aabb<T>(const triangle<T,2>& triangle, T& x1, T& y1, T& x2, T& y2);\
      template inline void aabb<T>(const rectangle<T>& rectangle, T& x1, T& y1, T& x2, T& y2);\
      template inline void aabb<T>(const quadix<T,2>& quadix, T& x1, T& y1, T& x2, T& y2);\
      template inline void aabb<T>(const circle<T>& circle, T& x1, T& y1, T& x2, T& y2);\
      template inline void aabb<T>(const polygon<T,2>& polygon,T& x1, T& y1, T& x2, T& y2);\
      template inline box<T,3> aabb<T>(const segment<T,3>& segment);\
      template inline box<T,3> aabb<T>(const triangle<T,3>& triangle);\
      template inline box<T,3> aabb<T>(const box<T,3>& rectangle);\
      template inline box<T,3> aabb<T>(const quadix<T,3>& quadix);\
      template inline box<T,3> aabb<T>(const sphere<T>& sphere);\
      template inline box<T,3> aabb<T>(const polygon<T,3>& polygon);\
      template inline void aabb<T>(const segment<T,3>& segment,T& x1, T& y1, T& z1, T& x2, T& y2, T& z2);\
      template inline void aabb<T>(const triangle<T,3>& triangle, T& x1, T& y1, T& z1, T& x2, T& y2, T& z2);\
      template inline void aabb<T>(const box<T,3>& box, T& x1, T& y1, T& z1, T& x2, T& y2, T& z2);\
      template inline void aabb<T>(const quadix<T,3>& quadix, T& x1, T& y1, T& z1, T& x2, T& y2, T& z2);\
      template inline void aabb<T>(const sphere<T>& sphere, T& x1, T& y1, T& z1, T& x2, T& y2, T& z2);\
      template inline void aabb<T>(const polygon<T,3>& polygon,T& x1, T& y1, T& z1, T& x2, T& y2, T& z2);\
      template inline rectangle<T> update_rectangle<T>(const rectangle<T>& rectangle, point2d<T>& point);\
      template inline box<T,3> update_box<T>(const box<T,3>& box, point3d<T>& point);\
      template inline circle<T> update_circle<T>(const circle<T>& circle, point2d<T>& point);\
      template inline sphere<T> update_sphere<T>(const sphere<T>& sphere, point3d<T>& point);\
      template inline point2d<T> generate_point_on_segment<T>(const segment<T,2>& segment, const T& t);\
      template inline point3d<T> generate_point_on_segment<T>(const segment<T,3>& segment, const T& t);\
      template inline point2d<T> generate_point_on_ray<T>(const ray<T,2>& ray, const T& t);\
      template inline point3d<T> generate_point_on_ray<T>(const ray<T,3>& ray, const T& t);\
      template inline T generate_random_value<T>(const T& range);\
      template inline point2d<T> generate_random_point<T>(const T& dx, const T& dy);\
      template inline point3d<T> generate_random_point<T>(const T& dx, const T& dy, const T& dz);\
      template inline point2d<T> generate_random_point<T>(const segment<T,2>& segment);\
      template inline point3d<T> generate_random_point<T>(const segment<T,3>& segment);\
      template inline point2d<T> generate_random_point<T>(const triangle<T,2>& triangle);\
      template inline point3d<T> generate_random_point<T>(const triangle<T,3>& triangle);\
      template inline point2d<T> generate_random_point<T>(const quadix<T,2>& quadix);\
      template inline point3d<T> generate_random_point<T>(const quadix<T,3>& quadix);\
      template inline point2d<T> generate_random_point<T>(const rectangle<T>& rectangle);\
      template inline point3d<T> generate_random_point<T>(const box<T,3>& box);\
      template inline void generate_random_points<T>(const T& x1, const T& y1, const T& x2, const T& y2, const std::size_t& point_count, OutputIterator2d out);\
      template inline void generate_random_points<T>(const rectangle<T>& rectangle, const std::size_t& point_count, OutputIterator2d out);\
      template inline void generate_random_points<T>(const segment<T,2>& segment, const std::size_t& point_count, OutputIterator2d out);\
      template inline void generate_random_points<T>(const segment<T,3>& segment, const std::size_t& point_count, OutputIterator2d out);\
      template inline void generate_random_points<T>(const triangle<T,2>& triangle, const std::size_t& point_count, OutputIterator2d out);\
      template inline void generate_random_points<T>(const triangle<T,3>& triangle, const std::size_t& point_count, OutputIterator2d out);\
      template inline void generate_random_points<T>(const quadix<T,2>& quadix, const std::size_t& point_count, OutputIterator2d out);\
      template inline void generate_random_points<T>(const quadix<T,3>& quadix, const std::size_t& point_count, OutputIterator2d out);\
      template inline void generate_random_points<T>(const circle<T>& circle, const std::size_t& point_count, OutputIterator2d out);\
      template inline void generate_random_object<T>(const T& x1, const  T& y1, const  T& x2, const  T& y2, segment<T,2>& segment);\
      template inline void generate_random_object<T>(const T& x1, const T& y1, const T& x2, const T& y2, rectangle<T>& rectangle);\
      template inline void generate_random_object<T>(const T& x1, const T& y1, const T& x2, const T& y2, triangle<T,2>& triangle);\
      template inline void generate_random_object<T>(const T& x1, const T& y1, const T& x2, const T& y2, quadix<T,2>& quadix);\
      template inline void generate_random_object<T>(const T& x1, const T& y1, const T& x2, const T& y2, circle<T>& circle);\
      template inline void generate_random_object<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, box<T,3>& box);\
      template inline triangle<T,2> right_shift<T>(const triangle<T,2>& triangle, const std::size_t& shift);\
      template inline triangle<T,3> right_shift<T>(const triangle<T,3>& triangle, const std::size_t& shift);\
      template inline quadix<T,2> right_shift<T>(const quadix<T,2>& quadix, const std::size_t& shift);\
      template inline quadix<T,3> right_shift<T>(const quadix<T,3>& quadix, const std::size_t& shift);\
      template inline T vector_norm<T>(const vector2d<T>& v);\
      template inline T vector_norm<T>(const vector3d<T>& v);\
      template inline vector2d<T> normalize<T>(const vector2d<T>& v);\
      template inline vector3d<T> normalize<T>(const vector3d<T>& v);\
      template inline vector2d<T> perpendicular<T>(const vector2d<T>& v);\
      template inline vector3d<T> perpendicular<T>(const vector3d<T>& v);\
      template inline vector2d<T> operator+<T>(const vector2d<T>& v1, const vector2d<T>& v2);\
      template inline vector3d<T> operator+<T>(const vector3d<T>& v1, const vector3d<T>& v2);\
      template inline vector2d<T> operator-<T>(const vector2d<T>& v1, const vector2d<T>& v2);\
      template inline vector3d<T> operator-<T>(const vector3d<T>& v1, const vector3d<T>& v2);\
      template inline T operator*<T>(const vector2d<T>& v1, const vector2d<T>& v2);\
      template inline vector3d<T> operator*<T>(const vector3d<T>& v1, const vector3d<T>& v2);\
      template inline T dot_product<T>(const vector2d<T>& v1, const vector2d<T>& v2);\
      template inline T dot_product<T>(const vector3d<T>& v1, const vector3d<T>& v2);\
      template inline T perpendicular_product<T>(const vector2d<T>& v1, const vector2d<T>& v2);\
      template inline T triple_product<T>(const vector3d<T>& v1, const vector3d<T>& v2, const vector3d<T>& v3);\
      template inline vector2d<T> operator*<T>(const vector2d<T>& v1, const T& scale);\
      template inline vector3d<T> operator*<T>(const vector3d<T>& v1, const T& scale);\
      template inline vector2d<T> operator*<T>(const T& scale, const vector2d<T>& v1);\
      template inline vector3d<T> operator*<T>(const T& scale, const vector3d<T>& v1);\
      template inline point2d<T> operator*<T>(const point2d<T>& point, const T& scale);\
      template inline point3d<T> operator*<T>(const point3d<T>& point, const T& scale);\
      template inline point2d<T> operator*<T>(const T& scale, const point2d<T>& point);\
      template inline point3d<T> operator*<T>(const T& scale, const point3d<T>& point);\
      template inline point2d<T> operator+<T>(const point2d<T>& point, const vector2d<T>& v);\
      template inline point2d<T> operator+<T>(const vector2d<T>& v, const point2d<T>& point);\
      template inline point3d<T> operator+<T>(const point3d<T>& point, const vector3d<T>& v);\
      template inline point3d<T> operator+<T>(const vector3d<T>& v, const point3d<T>& point);\
      template inline vector2d<T> operator-<T>(const point2d<T>& p1, const point2d<T>& p2);\
      template inline vector3d<T> operator-<T>(const point3d<T>& p1, const point3d<T>& p2);\
      template inline point2d<T> operator+<T>(const point2d<T>& p1, const point2d<T>& p2);\
      template inline point3d<T> operator+<T>(const point3d<T>& p1, const point3d<T>& p2);\
      template inline bool is_equal<T>(const T& val1, const T& val2, const T& epsilon);\
      template inline bool is_equal<T>(const point2d<T>& point1, const point2d<T>& point2, const T& epsilon);\
      template inline bool is_equal<T>(const point3d<T>& point1, const point3d<T>& point2, const T& epsilon);\
      template inline bool is_equal<T>(const T& val1, const T& val2);\
      template inline bool is_equal<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline bool is_equal<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline bool is_equal<T>(const rectangle<T>& rectangle1, const rectangle<T>& rectangle2);\
      template inline bool is_equal<T>(const circle<T>& circle1, const circle<T>& circle2);\
      template inline bool is_equal<T>(const box<T,3>& box1, const box<T,3>& box2);\
      template inline bool is_equal<T>(const sphere<T>& sphere1, const sphere<T>& sphere2);\
      template inline bool not_equal<T>(const T& val1, const T& val2, const T& epsilon);\
      template inline bool not_equal<T>(const point2d<T>& point1, const point2d<T>& point2, const T& epsilon);\
      template inline bool not_equal<T>(const point3d<T>& point1, const point3d<T>& point2, const T& epsilon);\
      template inline bool not_equal<T>(const T& val1, const T& val2);\
      template inline bool not_equal<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline bool not_equal<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline bool less_than_or_equal<T>(const T& val1, const T& val2, const T& epsilon);\
      template inline bool less_than_or_equal<T>(const T& val1, const T& val2);\
      template inline bool greater_than_or_equal<T>(const T& val1, const T& val2, const T& epsilon);\
      template inline bool greater_than_or_equal<T>(const T& val1, const T& val2);\
      template inline bool operator< <T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline bool operator< <T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline bool operator><T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline bool operator><T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline bool operator==<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline bool operator==<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline bool is_degenerate<T>(const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline bool is_degenerate<T>(const segment<T,2>& segment);\
      template inline bool is_degenerate<T>(const line<T,2>& line);\
      template inline bool is_degenerate<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline bool is_degenerate<T>(const segment<T,3>& segment);\
      template inline bool is_degenerate<T>(const line<T,3>& line);\
      template inline bool is_degenerate<T>(const triangle<T,2>& triangle);\
      template inline bool is_degenerate<T>(const triangle<T,3>& triangle);\
      template inline bool is_degenerate<T>(const quadix<T,2>& quadix);\
      template inline bool is_degenerate<T>(const quadix<T,3>& quadix);\
      template inline bool is_degenerate<T>(const rectangle<T>& rectangle);\
      template inline bool is_degenerate<T>(const circle<T>& circle);\
      template inline bool is_degenerate<T>(const sphere<T>& sphere);\
      template inline bool is_degenerate<T>(const circular_arc<T>& arc);\
      template inline point2d<T> degenerate_point2d<T>();\
      template inline point3d<T> degenerate_point3d<T>();\
      template inline vector2d<T> degenerate_vector2d<T>();\
      template inline vector3d<T> degenerate_vector3d<T>();\
      template inline ray<T,2> degenerate_ray2d<T>();\
      template inline ray<T,3> degenerate_ray3d<T>();\
      template inline line<T,2> degenerate_line2d<T>();\
      template inline line<T,3> degenerate_line3d<T>();\
      template inline segment<T,2> degenerate_segment2d<T>();\
      template inline segment<T,3> degenerate_segment3d<T>();\
      template inline triangle<T,2> degenerate_triangle2d<T>();\
      template inline triangle<T,3> degenerate_triangle3d<T>();\
      template inline quadix<T,2> degenerate_quadix2d<T>();\
      template inline quadix<T,3> degenerate_quadix3d<T>();\
      template inline rectangle<T> degenerate_rectangle<T>();\
      template inline circle<T> degenerate_circle<T>();\
      template inline sphere<T> degenerate_sphere<T>();\
      template inline point2d<T> positive_infinite_point2d<T>();\
      template inline point2d<T> negative_infinite_point2d<T>();\
      template inline point3d<T> positive_infinite_point3d<T>();\
      template inline point3d<T> negative_infinite_point3d<T>();\
      template inline point2d<T> make_point<T>(const T& x, const T& y);\
      template inline point3d<T> make_point<T>(const T& x, const T& y, const T& z);\
      template inline point2d<T> make_point<T>(const point3d<T> point);\
      template inline point3d<T> make_point<T>(const point2d<T> point, const T& z);\
      template inline point2d<T> make_point<T>(const circle<T>& circle);\
      template inline point3d<T> make_point<T>(const sphere<T>& sphere);\
      template inline vector2d<T> make_vector<T>(const T& x, const T& y);\
      template inline vector3d<T> make_vector<T>(const T& x, const T& y, const T& z);\
      template inline vector2d<T> make_vector<T>(const vector3d<T> v);\
      template inline vector3d<T> make_vector<T>(const vector2d<T> v, const T& z);\
      template inline vector2d<T> make_vector<T>(const point2d<T> point);\
      template inline vector3d<T> make_vector<T>(const point3d<T> point);\
      template inline ray<T,2> make_ray<T>(const T& ox, const T& oy, const T& dir_x, const T& dir_y);\
      template inline ray<T,3> make_ray<T>(const T& ox, const T& oy, const T& oz, const T& dir_x, const T& dir_y, const T& dir_z);\
      template inline ray<T,2> make_ray<T>(const point2d<T>& origin, const vector2d<T>& direction);\
      template inline ray<T,3> make_ray<T>(const point3d<T>& origin, const vector3d<T>& direction);\
      template inline ray<T,2> make_ray<T>(const point2d<T>& origin, const T& bearing);\
      template inline curve_point<T,2> make_curve_point<T>(const T& x, const T& y, const T& t);\
      template inline curve_point<T,3> make_curve_point<T>(const T& x, const T& y, const T& z, const T& t);\
      template inline curve_point<T,2> make_curve_point<T>(const point2d<T>& point, const T& t);\
      template inline curve_point<T,3> make_curve_point<T>(const point3d<T>& point, const T& t);\
      template inline segment<T,2> make_segment<T>(const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline segment<T,3> make_segment<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline segment<T,2> make_segment<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline segment<T,3> make_segment<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline segment<T,2> make_segment<T>(const line<T,2>& line);\
      template inline segment<T,3> make_segment<T>(const line<T,3>& line);\
      template inline line<T,2> make_line<T>(const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline line<T,3> make_line<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline line<T,2> make_line<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline line<T,3> make_line<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline line<T,2> make_line<T>(const segment<T,2>& segment);\
      template inline line<T,3> make_line<T>(const segment<T,3>& segment);\
      template inline line<T,2> make_line<T>(const ray<T,2>& ray);\
      template inline line<T,3> make_line<T>(const ray<T,3>& ray);\
      template inline rectangle<T> make_rectangle<T>(const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline rectangle<T> make_rectangle<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline box<T,3> make_box<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2);\
      template inline box<T,3> make_box<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline triangle<T,2> make_triangle<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3);\
      template inline triangle<T,3> make_triangle<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3);\
      template inline triangle<T,2> make_triangle<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline triangle<T,3> make_triangle<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3);\
      template inline quadix<T,2> make_quadix<T>(const T& x1, const T& y1, const T& x2, const T& y2, const T& x3, const T& y3, const T& x4, const T& y4);\
      template inline quadix<T,3> make_quadix<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3, const T& x4, const T& y4, const T& z4);\
      template inline quadix<T,2> make_quadix<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3, const point2d<T>& point4);\
      template inline quadix<T,3> make_quadix<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3, const point3d<T>& point4);\
      template inline quadix<T,2> make_quadix<T>(const T& x1, const T& y1, const T& x2, const T& y2);\
      template inline quadix<T,2> make_quadix<T>(const rectangle<T>& rectangle);\
      template inline circle<T> make_circle<T>(const T& x, const T& y, const T& radius);\
      template inline circle<T> make_circle<T>(const point2d<T>& point, const T& radius);\
      template inline circle<T> make_circle<T>(const point2d<T>& point1, const point2d<T>& point2);\
      template inline circle<T> make_circle<T>(const point2d<T>& point1, const point2d<T>& point2, const point2d<T>& point3);\
      template inline circle<T> make_circle<T>(const triangle<T,2>& triangle);\
      template inline sphere<T> make_sphere<T>(const T& x, const T& y, const T& z, const T& radius);\
      template inline sphere<T> make_sphere<T>(const point3d<T>& point, const T& radius);\
      template inline sphere<T> make_sphere<T>(const point3d<T>& point1, const point3d<T>& point2);\
      template inline plane<T,3> make_plane<T>(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, const T& x3, const T& y3, const T& z3);\
      template inline plane<T,3> make_plane<T>(const point3d<T>& point1, const point3d<T>& point2, const point3d<T>& point3);\
      template inline plane<T,3> make_plane<T>(const triangle<T,3>& triangle);\
      template inline polygon<T,2> make_polygon<T>(const std::vector< point2d<T> >& point_list);\
      template inline polygon<T,3> make_polygon<T>(const std::vector< point3d<T> >& point_list);\
      template inline polygon<T,2> make_polygon<T>(const triangle<T,2>& triangle);\
      template inline polygon<T,2> make_polygon<T>(const quadix<T,2>& quadix);\
      template inline polygon<T,2> make_polygon<T>(const rectangle<T>& rectangle);\
      template inline polygon<T,2> make_polygon<T>(const circle<T>& circle, const unsigned int point_count);


   #define INSTANTIATE_WYKOBI_ND(T, D, OutputIterator)\
      template inline bool parallel<T,D>(const line<T,D>& line1, const line<T,D>& line2);\
      template inline bool parallel<T,D>(const segment<T,D>& segment1, const segment<T,D>& segment2);\
      template inline bool collinear<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2, const pointnd<T,D>& point3);\
      template inline bool robust_collinear<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2, const pointnd<T,D>& point3);\
      template inline bool is_point_collinear<T,D>(const segment<T,D>& segment, const pointnd<T,D>& point, const bool robust);\
      template inline bool perpendicular<T,D>(const line<T,D>& line1, const line<T,D>& line2);\
      template inline bool perpendicular<T,D>(const segment<T,D>& segment1, const segment<T,D>& segment2);\
      template inline bool intersect<T,D>(const segment<T,D>& segment1, const segment<T,D>& segment2, const T& fuzzy);\
      template inline bool intersect<T,D>(const line<T,D>& line1, const line<T,D>& line2, const T& fuzzy);\
      template inline pointnd<T,D> intersection_point<T,D>(const segment<T,D>& segment1, const segment<T,D>& segment2, const T& fuzzy);\
      template inline pointnd<T,D> intersection_point<T,D>(const line<T,D>& line1, const line<T,D>& line2, const T& fuzzy);\
      template inline T distance<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline T distance<T,D>(const pointnd<T,D>& point, const segment<T,D>& segment);\
      template inline T distance<T,D>(const pointnd<T,D>& point, const line<T,D>& line);\
      template inline T distance<T,D>(const segment<T,D>& segment1, const segment<T,D>& segment2);\
      template inline T distance<T,D>(const line<T,D>& line1, const line<T,D>& line2);\
      template inline T lay_distance<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline T lay_distance<T,D>(const pointnd<T,D>& point, const segment<T,D>& segment);\
      template inline T lay_distance<T,D>(const pointnd<T,D>& point, const line<T,D>& line);\
      template inline T lay_distance<T,D>(const segment<T,D>& segment1, const segment<T,D>& segment2);\
      template inline T lay_distance<T,D>(const line<T,D>& line1, const line<T,D>& line2);\
      template inline T manhattan_distance<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline T chebyshev_distance<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline T inverse_chebyshev_distance<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline bool point_in_box<T,D>(const pointnd<T,D>& point, const box<T,D>& box);\
      template inline bool point_in_sphere<T,D>(const pointnd<T,D>& point, const hypersphere<T,D>& sphere);\
      template inline pointnd<T,D> closest_point_on_segment_from_point<T,D>(const segment<T,D>& segment, const wykobi::segment<T,D>::PointType& point);\
      template inline pointnd<T,D> closest_point_on_line_from_point<T,D>(const line<T,D>& segment, const wykobi::line<T,D>::PointType& point);\
      template inline pointnd<T,D> closest_point_on_sphere_from_point<T,D>(const hypersphere<T,D>& sphere, const wykobi::hypersphere<T,D>::PointType& point);\
      template inline pointnd<T,D> closest_point_on_plane_from_point<T,D>(const plane<T,D>& plane, const wykobi::plane<T,D>::PointType& point);\
      template inline pointnd<T,D> closest_point_on_box_from_point<T,D>(const box<T,D>& box, const wykobi::box<T,D>::PointType& point);\
      template inline pointnd<T,D> project_point_t<T,D>(const pointnd<T,D>& source_point, const pointnd<T,D>& destination_point, const T& t);\
      template inline pointnd<T,D> project_point<T,D>(const pointnd<T,D>& source_point, const pointnd<T,D>& destination_point, const T& distance);\
      template inline pointnd<T,D> mirror<T,D>(const pointnd<T,D>& point, const line<T,D>& mirror_axis);\
      template inline segment<T,D> mirror<T,D>(const segment<T,D>& segment, const line<T,D>& mirror_axis);\
      template inline line<T,D> mirror<T,D>(const line<T,D>& line, const wykobi::line<T,D>& mirror_axis);\
      template inline box<T,D> mirror<T,D>(const box<T,D>& box, const line<T,D>& mirror_axis);\
      template inline triangle<T,D> mirror<T,D>(const triangle<T,D>& triangle, const line<T,D>& mirror_axis);\
      template inline quadix<T,D> mirror<T,D>(const quadix<T,D>& quadix, const line<T,D>& mirror_axis);\
      template inline hypersphere<T,D> mirror<T,D>(const hypersphere<T,D>& sphere, const line<T,D>& mirror_axis);\
      template inline polygon<T,D> mirror<T,D>(const polygon<T,D>& polygon, const line<T,D>& mirror_axis);\
      template inline segment<T,D> project_onto_axis<T,D>(const pointnd<T,D>& point, const line<T,D>& axis);\
      template inline segment<T,D> project_onto_axis<T,D>(const triangle<T,D>& triangle, const line<T,D>& axis);\
      template inline segment<T,D> project_onto_axis<T,D>(const box<T,D>& rectangle, const line<T,D>& axis);\
      template inline segment<T,D> project_onto_axis<T,D>(const quadix<T,D>& quadix, const line<T,D>& axis);\
      template inline segment<T,D> project_onto_axis<T,D>(const hypersphere<T,D>& sphere, const line<T,D>& axis);\
      template inline segment<T,D> project_onto_axis<T,D>(const polygon<T,D>& polygon, const line<T,D>& axis);\
      template inline T perimeter<T,D>(const triangle<T,D>& triangle);\
      template inline T perimeter<T,D>(const quadix<T,D>& quadix);\
      template inline T perimeter<T,D>(const polygon<T,D>& polygon);\
      template inline pointnd<T,D> generate_random_point<T,D>(const segment<T,D>& segment);\
      template inline pointnd<T,D> generate_random_point<T,D>(const triangle<T,D>& triangle);\
      template inline pointnd<T,D> generate_random_point<T,D>(const quadix<T,D>& quadix);\
      template inline pointnd<T,D> generate_random_point<T,D>(const box<T,D>& box);\
      template inline void generate_random_points<T,D,OutputIterator>(const box<T,D>& box, const std::size_t& point_count, OutputIterator out);\
      template inline void generate_random_points<T,D,OutputIterator>(const segment<T,D>& segment, const std::size_t& point_count, OutputIterator out);\
      template inline void generate_random_points<T,D,OutputIterator>(const triangle<T,D>& triangle, const std::size_t& point_count, OutputIterator out);\
      template inline void generate_random_points<T,D,OutputIterator>(const quadix<T,D>& quadix, const std::size_t& point_count, OutputIterator out);\
      template inline T vector_norm<T,D>(const vectornd<T,D>& v);\
      template inline vectornd<T,D> normalize<T,D>(const vectornd<T,D>& v);\
      template inline vectornd<T,D> operator+<T,D>(const vectornd<T,D>& v1, const vectornd<T,D>& v2);\
      template inline vectornd<T,D> operator-<T,D>(const vectornd<T,D>& v1, const vectornd<T,D>& v2);\
      template inline T dot_product<T,D>(const vectornd<T,D>& v1, const vectornd<T,D>& v2);\
      template inline vectornd<T,D> operator*<T,D>(const vectornd<T,D>& v1, const T& scale);\
      template inline vectornd<T,D> operator*<T,D>(const T& scale, const vectornd<T,D>& v1);\
      template inline pointnd<T,D> operator*<T,D>(const pointnd<T,D>& point, const T& scale);\
      template inline pointnd<T,D> operator*<T,D>(const T& scale, const pointnd<T,D>& point);\
      template inline pointnd<T,D> operator+<T,D>(const pointnd<T,D>& point, const vectornd<T,D>& v);\
      template inline pointnd<T,D> operator+<T,D>(const vectornd<T,D>& v, const pointnd<T,D>& point);\
      template inline vectornd<T,D> operator-<T,D>(const pointnd<T,D>& p1, const pointnd<T,D>& p2);\
      template inline pointnd<T,D> operator+<T,D>(const pointnd<T,D>& p1, const pointnd<T,D>& p2);\
      template inline bool operator < <T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline bool operator > <T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline bool operator==<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline bool is_equal<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2, const T& epsilon);\
      template inline bool is_equal<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline bool not_equal<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2, const T& epsilon);\
      template inline bool not_equal<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline pointnd<T,D> degenerate_pointnd<T,D>();\
      template inline vectornd<T,D> degenerate_vectornd<T,D>();\
      template inline ray<T,D> degenerate_raynd<T,D>();\
      template inline line<T,D> degenerate_linend<T,D>();\
      template inline segment<T,D> degenerate_segmentnd<T,D>();\
      template inline triangle<T,D> degenerate_trianglend<T,D>();\
      template inline quadix<T,D> degenerate_quadixnd<T,D>();\
      template inline box<T,D> degenerate_box<T,D>();\
      template inline hypersphere<T,D> degenerate_hypersphere<T,D>();\
      template inline pointnd<T,D> positive_infinite_pointnd<T,D>();\
      template inline pointnd<T,D> negative_infinite_pointnd<T,D>();\
      template inline void swap<T,D>(pointnd<T,D>& point1, pointnd<T,D>& point2);\
      template inline vectornd<T,D> make_vector<T,D>(const pointnd<T,D>& point);\
      template inline ray<T,D> make_ray<T,D>(const pointnd<T,D>& origin, const vectornd<T,D>& direction);\
      template inline segment<T,D> make_segment<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline line<T,D> make_line<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline box<T,D> make_box<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline triangle<T,D> make_triangle<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2, const pointnd<T,D>& point3);\
      template inline quadix<T,D> make_quadix<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2, const pointnd<T,D>& point3, const pointnd<T,D>& point4);\
      template inline hypersphere<T,D> make_sphere<T,D>(const pointnd<T,D>& point, const T& radius);\
      template inline hypersphere<T,D> make_sphere<T,D>(const pointnd<T,D>& point1, const pointnd<T,D>& point2);\
      template inline polygon<T,D> make_polygon<T,D>(const std::vector< pointnd<T,D> >& point_list);\
      template inline polygon<T,D> make_polygon<T,D>(const triangle<T,D>& triangle);\
      template inline polygon<T,D> make_polygon<T,D>(const quadix<T,D>& quadix);

   #define INSTANTIATE_WYKOBI_MATH(T)\
      template inline T sqr<T>(const T& val);\
      template inline T abs<T>(const T& value);\
      template inline T max<T>(const T& value1, const T& value2);\
      template inline T min<T>(const T& value1, const T& value2);\
      template inline T infinity<T>();\
      template inline T sin<T>(const T& value);\
      template inline T cos<T>(const T& value);\
      template inline T tan<T>(const T& value);\
      template inline T atan<T>(const T& value);\
      template inline T approx_sin<T>(T angle);\
      template inline T approx_cos<T>(T angle);\
      template inline T approx_tan<T>(T angle);\
      template inline T clamp(const T& value, const T& low, const T& high);


   #define INSTANTIATE_WYKOBI_UTILITIES_1(T)\
      template std::ostream& operator<< <T>(std::ostream& os, const point2d<T>& point);\
      template std::ostream& operator<< <T>(std::ostream& os, const point3d<T>& point);\
      template std::ostream& operator<< <T>(std::ostream& os, const ray<T,2>& ray);\
      template std::ostream& operator<< <T>(std::ostream& os, const ray<T,3>& ray);\
      template std::ostream& operator<< <T>(std::ostream& os, const vector2d<T>& v);\
      template std::ostream& operator<< <T>(std::ostream& os, const vector3d<T>& v);\
      template std::ostream& operator<< <T>(std::ostream& os, const circle<T>& circle);\
      template std::ostream& operator<< <T>(std::ostream& os, const sphere<T>& sphere);\
      template std::ostream& operator<< <T>(std::ostream& os, const rectangle<T>& rectangle);\
      template std::ostream& operator<< <T>(std::ostream& os, const box<T,3>& box);


   #define INSTANTIATE_WYKOBI_UTILITIES_2(T,D)\
      template std::ostream& operator<< <T,D>(std::ostream& os, const segment<T,D>& segment);\
      template std::ostream& operator<< <T,D>(std::ostream& os, const line<T,D>& line);\
      template std::ostream& operator<< <T,D>(std::ostream& os, const triangle<T,D>& triangle);\
      template std::ostream& operator<< <T,D>(std::ostream& os, const quadix<T,D>& quadix);

   typedef wykobi::point2d<float>*  flt_pnt_2d;
   typedef wykobi::point2d<double>* dbl_pnt_2d;

   typedef wykobi::point3d<float>*  flt_pnt_3d;
   typedef wykobi::point3d<double>* dbl_pnt_3d;

   INSTANTIATE_WYKOBI(float,flt_pnt_2d,flt_pnt_3d,flt_pnt_2d,flt_pnt_3d)
   INSTANTIATE_WYKOBI(double,dbl_pnt_2d,dbl_pnt_3d,dbl_pnt_2d,dbl_pnt_3d)


   /*
   typedef pointnd<float, 4> pointnd_flt_4;
   typedef pointnd<float, 5> pointnd_flt_5;
   typedef pointnd<float, 6> pointnd_flt_6;
   typedef pointnd<float, 7> pointnd_flt_7;
   typedef pointnd<float, 8> pointnd_flt_8;
   typedef pointnd<float, 9> pointnd_flt_9;
   typedef pointnd<float,10> pointnd_flt_10;

   typedef pointnd<double, 4> pointnd_dbl_4;
   typedef pointnd<double, 5> pointnd_dbl_5;
   typedef pointnd<double, 6> pointnd_dbl_6;
   typedef pointnd<double, 7> pointnd_dbl_7;
   typedef pointnd<double, 8> pointnd_dbl_8;
   typedef pointnd<double, 9> pointnd_dbl_9;
   typedef pointnd<double,10> pointnd_dbl_10;

   INSTANTIATE_WYKOBI_ND(float, 4,  pointnd_flt_4)
   INSTANTIATE_WYKOBI_ND(float, 5,  pointnd_flt_5)
   INSTANTIATE_WYKOBI_ND(float, 6,  pointnd_flt_6)
   INSTANTIATE_WYKOBI_ND(float, 7,  pointnd_flt_7)
   INSTANTIATE_WYKOBI_ND(float, 8,  pointnd_flt_8)
   INSTANTIATE_WYKOBI_ND(float, 9,  pointnd_flt_9)
   INSTANTIATE_WYKOBI_ND(float,10, pointnd_flt_10)

   INSTANTIATE_WYKOBI_ND(double, 4,  pointnd_dbl_4)
   INSTANTIATE_WYKOBI_ND(double, 5,  pointnd_dbl_5)
   INSTANTIATE_WYKOBI_ND(double, 6,  pointnd_dbl_6)
   INSTANTIATE_WYKOBI_ND(double, 7,  pointnd_dbl_7)
   INSTANTIATE_WYKOBI_ND(double, 8,  pointnd_dbl_8)
   INSTANTIATE_WYKOBI_ND(double, 9,  pointnd_dbl_9)
   INSTANTIATE_WYKOBI_ND(double,10, pointnd_dbl_10)
   */

   INSTANTIATE_WYKOBI_MATH(float)
   INSTANTIATE_WYKOBI_MATH(double)

   INSTANTIATE_WYKOBI_UTILITIES_1(float)
   INSTANTIATE_WYKOBI_UTILITIES_1(double)

   INSTANTIATE_WYKOBI_UTILITIES_2(float,2)
   INSTANTIATE_WYKOBI_UTILITIES_2(double,2)

   INSTANTIATE_WYKOBI_UTILITIES_2(float,3)
   INSTANTIATE_WYKOBI_UTILITIES_2(double,3)

} // wykobi namespace


#endif
