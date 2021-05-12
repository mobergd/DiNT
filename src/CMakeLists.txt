if(USE_TINKER) # NEEDS FIX
    set(POT_DIR "${CMAKE_SOURCE_DIR}/dd_pot")
else()
    set(POT_DIR "${CMAKE_SOURCE_DIR}/pot")
endif(USE_TINKER)

set(src_SOURCES 
 addrot.F 	airy.f		andersen.F 	angmom.f 	ange.f 		atomdiatom.F 	bootstrap.f 
 bsstep.f 	carttojac.f 	checkfrus.f 	checkhop.F 	cn.f 		decocheck.F
 derivs.f 	detwell.f 	diamin.f 	diamin2.f 	diapot.f 	diapot2.f 
 diatom.F 	dint.F 		driver.f 	ejscale.f 	elecdeco.f 	endmpi.f	ems.F
 ewkb.f 	ewkb2.f		finalstate.f 	getdvec2.f 	getgrad.f 	getpem.f 
 getrel.f 	getrelbz2.f	getrho.f 	gettemp.f 	header.f 	hessems.F 	hop.f 
 honey.f 	initelec.f 	initmol.f	initmpi.f	jointprint.F	lindemann.f	mmid.f 
 noang.f 	normod.f 	nmturn.f 	nmpot.f 	orient.F 
 period.f 	period2.f 	pjsplit.f 	getplz.f 	popnorm.F 	preatomdiatom.f 
 prediatom.f 	premol.f 	preptraj.f 	pzextr.f 	radialdist.f 	ranclu.F 
 rancol.F 	rangaussian.F 	ranj.F 		ranno.F 	ranrot.F 	ransamp.F 
 ranstate.F 	rantherm.F 	rarray.f 	readin.f 	rk4.f 
 rtox.f 	spin.F 		stodeco.F 	takestep.f 	timing.f 	torphi.f 
 turn.f 	turn2.f 	vibwells.F 	vwkb.f 		vwkb2.f 	xptoy.f
 ${POT_DIR}/${POT}.f
)

if(USE_MPI)
  list(APPEND EXTRA_INCLUDES "${MPI_INCLUDE_PATH}")
  list(APPEND EXTRA_LIBS "${MPI_Fortran_LIBRARIES}")
endif()

#list(APPEND EXTRA_LIBS "${SPRNG_DIR}/lib")
list(APPEND EXTRA_LIBS "${SPRNG_DIR}/lib/libcmrg.a" "${SPRNG_DIR}/lib/liblcg64.a" "${SPRNG_DIR}/lib/liblcg.a" "${SPRNG_DIR}/lib/liblfg.a" "${SPRNG_DIR}/lib/libsprngtest.a")

add_library(objs OBJECT ${src_SOURCES})
list(APPEND EXTRA_INCLUDES "${SPRNG_DIR}/SRC")
list(APPEND EXTRA_INCLUDES "${CMAKE_SOURCE_DIR}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
#target_link_libraries(objs PRIVATE ${EXTRA_LIBS})
#target_include_directories(objs PRIVATE ${EXTRA_INCLUDES})
target_link_libraries(objs PUBLIC ${EXTRA_LIBS})
target_include_directories(objs PUBLIC ${EXTRA_INCLUDES})
