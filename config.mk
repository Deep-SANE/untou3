ifdef UNTOU3_DIR
CXXFLAGS += -I$(UNTOU3_DIR)/include
endif

ifdef HAVE_BOOST
CXXFLAGS += -DHAVE_BOOST
ifdef BOOST_ROOT
CXXFLAGS += -I$(BOOST_ROOT)/include
endif
endif
