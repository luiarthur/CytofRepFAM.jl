all: build

build:
	latexmk -pdf

clean:
	latexmk -C

watch:
	watchmedo shell-command \
			--patterns="*.tex;" \
			--recursive \
			--command='latexmk -pdf && latexmk -c' \
			--wait \
			.
