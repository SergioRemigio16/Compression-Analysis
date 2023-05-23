#!/bin/bash

# Install Node.js for coc.nvim
curl -sL https://deb.nodesource.com/setup_14.x | bash - && \
    apt-get install -y nodejs

# Install Vundle for Vim
git clone https://github.com/VundleVim/Vundle.vim.git ~/.vim/bundle/Vundle.vim

# Convert .vimrc to Unix format and install Vim plugins
dos2unix /root/.vimrc && vim +PluginInstall +qall

# Compile coc.nvim
cd ~/.vim/bundle/coc.nvim && yarn install && yarn build