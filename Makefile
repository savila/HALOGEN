EXEC   = 2LPT-HALOGEN

$(EXEC):
		$(MAKE) -C src $(EXEC);\
		mv -f src/$(EXEC) $(EXEC)

2LPT:
		$(MAKE) -C src 2LPT;\
		mv -f src/2LPT 2LPT


halogen:
		$(MAKE) -C src halogen;\
		mv -f src/halogen halogen

fit:
		$(MAKE) -C src fit;\
		mv -f src/fit fit

clean:
		$(MAKE) -C src clean;\
		rm -rf halogen
		rm -f fit
		rm -f 2LPT
		rm -f $(EXEC)
