
set encoding=utf-8

" This is a sample .vimrc file for C++ development
" You will need to install some plugins to fully utilize this

set nocompatible              " be iMproved, required
filetype off                  " required

" set the runtime path to include Vundle and initialize
set rtp+=~/.vim/bundle/Vundle.vim
call vundle#begin()

" let Vundle manage Vundle, required
Plugin 'VundleVim/Vundle.vim'

" Add all your plugins here (note older versions of Vundle used Bundle instead of Plugin)
Plugin 'preservim/nerdtree'
Plugin 'vim-syntastic/syntastic'
Plugin 'nvie/vim-flake8'
Plugin 'neoclide/coc.nvim', {'branch': 'release'}
Plugin 'jiangmiao/auto-pairs'
Plugin 'preservim/nerdcommenter'
Plugin 'plasticboy/vim-markdown'
Plugin 'airblade/vim-gitgutter'
Plugin 'jeffkreeftmeijer/vim-hlslens'
Plugin 'altercation/vim-colors-solarized'  " NEW: Solarized color scheme

" All of your Plugins must be added before the following line
call vundle#end()            " required
filetype plugin indent on    " required

" General configurations
syntax on                    " syntax highlighting
set number                   " line numbers
set autoindent               " auto indentation
set smartindent              " smart indentation
set tabstop=4                " a tab is four spaces
set shiftwidth=4             " number of spaces to use for autoindenting
set expandtab                " use spaces instead of tabs
set cursorline

" NERDTree configurations
map <C-n> :NERDTreeToggle<CR>
let NERDTreeIgnore=['\~$', '\.pyc$', '\.swp$']

" Syntastic configurations
set statusline+=%#warningmsg#
set statusline+=%{SyntasticStatuslineFlag()}
set statusline+=%*

let g:syntastic_always_populate_loc_list = 1
let g:syntastic_auto_loc_list = 1
let g:syntastic_check_on_open = 1
let g:syntastic_check_on_wq = 0

" Auto Pairs configurations
let g:AutoPairsFlyMode = 1

" NERD Commenter configurations
let g:NERDSpaceDelims = 1

" GitGutter configurations
let g:gitgutter_enabled = 1
let g:gitgutter_eager = 0

" Highlight character under cursor
" HighlightSearchLens will only highlight the current match
let g:hlslens#enable = 1

" colorscheme elflord

" End of .vimrc file
