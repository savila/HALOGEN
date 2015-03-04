EXEC   = 2LPT-HALOGEN

$(EXEC):
		$(MAKE) -C src $(EXEC);\
		mv -f src/$(EXEC) $(EXEC)


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
		rm -f $(EXEC)
