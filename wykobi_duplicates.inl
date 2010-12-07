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


#include "wykobi.hpp"
#include "wykobi_algorithm.hpp"


namespace wykobi
{
   namespace algorithm
   {

      template<typename T>
      struct remove_duplicates
      {
      public:
         template<typename InputIterator, typename OutputIterator>
         remove_duplicates(InputIterator begin, InputIterator end, OutputIterator out)
         {
            std::sort(begin,end);
            T previous = (*begin);
            (*out++) = (*begin);
            for(InputIterator it = (begin + 1); it != end; ++it)
            {
               if ((*it) > previous)
               {
                  (*out++) = (*it);
                  previous = (*it);
               }
            }
         }
      };

   } // namespace wykobi::algorithm

} // namespace wykobi
