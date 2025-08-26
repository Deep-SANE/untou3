ifdef UNTOU3_DIR
CXXFLAGS += -I$(UNTOU3_DIR)/include
endif

ifdef U3LIB_DIR
CXXFLAGS += -I$(U3LIB_DIR)/include
endif

ifdef U4LIB_DIR
CXXFLAGS += -I$(U4LIB_DIR)/include
endif

ifdef QCX_UTILS_DIR
CXXFLAGS += -I$(QCX_UTILS_DIR)/include
endif

ifdef HAVE_BOOST
CXXFLAGS += -DHAVE_BOOST
ifdef BOOST_ROOT
CXXFLAGS += -I$(BOOST_ROOT)/include
endif
endif

