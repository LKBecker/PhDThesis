   M S B A M L   `   `   `  ���� G  CFastQCRipper, Version=1.0.0.0, Culture=neutral, PublicKeyToken=null+clr-namespace:FastQCRipperFastQCRipper      FastQCRipper.MainWindow   W SPresentationCore, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35R NWindowsBase, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35\ XPresentationFramework, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35D 9http://schemas.microsoft.com/winfx/2006/xaml/presentation    R NSystem.Xaml, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml     @d2http://schemas.microsoft.com/expression/blend/2008     Jmc;http://schemas.openxmlformats.org/markup-compatibility/2006     ,localclr-namespace:FastQCRipper       � Title$  Fast QC Ripper�� � SizeToContent$ WidthAndHeight=� � 
ResizeMode$ CanResizeWithGrip=�$��Center=�$��Center=�.��� -    �Name$ MainGrid��$	��521��$	��105�� J� RenderTransformOrigin$ 0.5,0.5&�$	��Top=�.+��� -    ��Name$! ButtonSelectSourceDirectory��$��	Select...��$
��Left=�$��	436,9,0,0q�$	��Top=�$��75���� -    ��Name$ TextBoxSourceDirectory��$
��Left=�$��23��$��	10,10,0,0q� �� TextWrapping$
 Wrap=�$��(No directory selected)��$	��Top=�$	��400��4  .System.Windows.Automation.AutomationProperties   Name$ TextBox.SourceDirectory��	   IsRowHeader	 . 
 �� 
IsReadOnly
 .  �� MaxLines$ 1�� �� MinLines$ 1��
���   ��� -   $&  ButtonSelectDestinationDirectory��$��	Select...��$
��Left=�$��
436,37,0,0q�$	��Top=�$��75���� -   $! TextBoxDestinationDirectory��$
��Left=�$��23��$��	10,38,0,0q�$
 Wrap=�$��(No directory selected)��$	��Top=�$	��400��$ 1��
 . 
���   ��� -   $ ButtomCompute��$��Extract��$
��Left=�$��
436,75,0,0q�$	��Top=�$��75��� -    �Name$ ProgressBarZIPRip��$��Center=�$��19��$��0q�$	��400��$��Center=�$y�Hidden=��� -    ��Name$ LabelProgress��$'��!Processing File X of Y (ZZZ%) ...��$��Center=�$��0,65,0,0q�$	��Top=�$y�Hidden=�9� $��Options��$��61��$��10,75,111,-31q�$��Bottom=� 9� 
IsExpanded .  .���� $
��Left=�$��51��$��10,0,0,0q�$	��380�� �� Orientation$ 
Horizontal=�.+��� $��File search pattern��$��23��$��0q�$	��110��$
��Left=�$��Center=��� $
��Left=�$��23��$
 Wrap=�$��TextBox��$��Center=�$	��240��$��
10,14,0,14q�