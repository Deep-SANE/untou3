include $(QCX_CONFIG_DIR)/make.d/*.mk

export 

.PHONY: all
all: 
	$(MAKE) -C src

.PHONY: clean
clean:
	$(MAKE) -C src clean
